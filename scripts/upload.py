#!/usr/bin/env python3
"""一键编译并上传 pywfn 到 PyPI。

流程：
  1. 版本冲突检查：读取本地版本号，查询 PyPI 上已有的版本。
     - 若本地版本已存在于 PyPI，或本地版本不高于线上最新版本，则判定为冲突，
       提示先升级版本号，并直接退出（不会编译、不会上传）。
  2. 同时编译 Windows 与 Linux 两个平台的 wheel：
     - Windows：调用本机（.venv 里）的 maturin 直接编译。
     - Linux  ：通过 WSL 调用其中的 maturin 编译。
  3. 使用 twine 将 dist/ 下生成的 wheel 上传到 PyPI。

示例（请在 Windows 侧用项目自带的 venv 运行）：
  .venv/Scripts/python scripts/upload.py               # 检查 -> 编译 -> 二次确认后上传
  .venv/Scripts/python scripts/upload.py --yes         # 上传前不再确认
  .venv/Scripts/python scripts/upload.py --no-upload   # 只编译，不上传
  .venv/Scripts/python scripts/upload.py -r testpypi   # 上传到 TestPyPI
  .venv/Scripts/python scripts/upload.py --only linux  # 只编译 Linux

注意：本脚本需要在 Windows 侧运行（因为要通过 wsl.exe 调用 WSL 编译 Linux 版本）。
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import NoReturn

# 保证中文在 Windows 控制台也能正常输出
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]


# ------------------------------ 基本配置 ------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = PROJECT_ROOT / "pyproject.toml"
DIST_DIR = PROJECT_ROOT / "dist"

PYPI_JSON_URL = "https://pypi.org/pypi/{name}/json"
DEFAULT_WSL_DISTRO = "Ubuntu-22.04"

# WSL 内部的 cargo 目标目录：放到 Linux 家目录，避免和 Windows 的 target/ 相互污染
LINUX_TARGET_DIR = "$HOME/.cache/pywfn-linux-target"


def log(msg: str = "") -> None:
    print(msg, flush=True)


def die(msg: str, code: int = 1) -> NoReturn:
    log(f"错误: {msg}")
    sys.exit(code)


# ------------------------------ 元数据解析 ------------------------------
def read_project_meta() -> tuple[str, str]:
    """从 pyproject.toml 的 [project] 段解析出包名和版本号。"""
    if not PYPROJECT.exists():
        die(f"找不到 {PYPROJECT}")
    text = PYPROJECT.read_text(encoding="utf-8")
    m = re.search(r"(?ms)^\[project\]\s*$(.*?)(?=^\[|\Z)", text)
    block = m.group(1) if m else text
    name = re.search(r'(?m)^\s*name\s*=\s*"([^"]+)"', block)
    version = re.search(r'(?m)^\s*version\s*=\s*"([^"]+)"', block)
    if not name or not version:
        die("无法从 pyproject.toml 的 [project] 段解析出 name / version")
    return name.group(1), version.group(1)


def version_key(v: str) -> tuple[int, ...]:
    """把版本号转成可比较的元组，例如 '1.0.18' -> (1, 0, 18)。

    只取数字段，足够本项目这种简单版本号使用。
    """
    parts = re.findall(r"\d+", v)
    return tuple(int(x) for x in parts) if parts else (0,)


# ------------------------------ 版本冲突检查 ------------------------------
def fetch_pypi_releases(name: str):
    """获取 PyPI 上已发布的版本字典；若包不存在则返回 None。"""
    url = PYPI_JSON_URL.format(name=name)
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "pywfn-release-script"},
    )
    try:
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.load(resp)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        die(f"查询 PyPI 失败: HTTP {e.code} {e.reason}")
    except urllib.error.URLError as e:
        die(f"查询 PyPI 失败（网络问题）: {e.reason}")
    return data.get("releases", {})


def check_version_conflict(name: str, version: str) -> None:
    """检查本地版本与 PyPI 已有版本是否冲突，冲突则提示升级并退出。"""
    releases = fetch_pypi_releases(name)
    if not releases:
        log(f"[版本检查] PyPI 上还没有 {name} 的任何版本，跳过冲突检查。")
        return

    versions = list(releases.keys())
    latest = max(versions, key=version_key)
    log(f"[版本检查] 本地版本: {version} | PyPI 最新版本: {latest}（共 {len(versions)} 个已发布版本）")

    if version in releases:
        die(
            f"版本 {version} 已经存在于 PyPI，不能重复上传。\n"
            f"      请先升级版本号（当前最新为 {latest}），修改 pyproject.toml 的 version 后再运行本脚本。",
            code=2,
        )

    if version_key(version) <= version_key(latest):
        die(
            f"本地版本 {version} 不高于 PyPI 上的最新版本 {latest}，直接上传会造成版本混乱。\n"
            f"      请先把版本号升级到大于 {latest} 的值，再运行本脚本。",
            code=2,
        )

    log("[版本检查] 通过，可以继续编译。")


# ------------------------------ 工具定位 ------------------------------
def find_exe(name: str) -> str:
    """在 venv / PATH 中查找可执行文件。"""
    candidates = [
        Path(sys.executable).parent / f"{name}.exe",  # 当前 python 所在目录（通常是 .venv/Scripts）
        PROJECT_ROOT / ".venv" / "Scripts" / f"{name}.exe",  # 项目自带的 venv
    ]
    for cand in candidates:
        if cand.exists():
            return str(cand)
    found = shutil.which(name)
    if found:
        return found
    die(f"找不到可执行文件 {name}（请确认已安装，或位于项目 .venv 中）")


def to_wsl_path(path: Path) -> str:
    """把 Windows 路径 (E:\\code\\pywfn) 转成 WSL 路径 (/mnt/e/code/pywfn)。"""
    s = str(Path(path).resolve())
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", s)
    if not m:
        return s.replace("\\", "/")
    drive, rest = m.group(1).lower(), m.group(2).replace("\\", "/")
    return f"/mnt/{drive}/{rest}"


def run(cmd, cwd: Path | None = None) -> None:
    log("")
    log("$ " + " ".join(str(c) for c in cmd))
    # 强制子进程（特别是 twine）使用 UTF-8，避免在 GBK 控制台下
    # rich 进度条渲染中文/符号时抛出 UnicodeEncodeError。
    env = dict(os.environ)
    env.setdefault("PYTHONUTF8", "1")
    env.setdefault("PYTHONIOENCODING", "utf-8")
    result = subprocess.run(
        [str(c) for c in cmd], cwd=str(cwd) if cwd else None, check=False, env=env
    )
    if result.returncode != 0:
        die(f"命令执行失败（退出码 {result.returncode}）")


# ------------------------------ 编译 ------------------------------
def clean_dist() -> None:
    """清空 dist/ 下的旧 wheel，避免把过期的包一起上传。"""
    DIST_DIR.mkdir(parents=True, exist_ok=True)
    removed = 0
    for f in DIST_DIR.glob("*.whl"):
        f.unlink()
        removed += 1
    if removed:
        log(f"[清理] 已删除 dist/ 下 {removed} 个旧 wheel")


def build_windows(out_dir: Path) -> None:
    log("== 编译 Windows wheel ==")
    maturin = find_exe("maturin")
    run([maturin, "build", "--release", "--out", str(out_dir)], cwd=PROJECT_ROOT)


def build_linux(out_dir: Path, distro: str) -> None:
    log("== 编译 Linux wheel（通过 WSL）==")
    wsl = shutil.which("wsl") or shutil.which("wsl.exe")
    if not wsl:
        die("找不到 wsl.exe，无法编译 Linux 版本（请在 Windows 侧运行本脚本）")
    proj = to_wsl_path(PROJECT_ROOT)
    out = to_wsl_path(out_dir)
    inner = (
        "set -e; "
        'export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"; '
        f'cd "{proj}"; '
        "maturin --version; "
        f'CARGO_TARGET_DIR="{LINUX_TARGET_DIR}" maturin build --release --out "{out}"'
    )
    run([wsl, "-d", distro, "--", "bash", "-lc", inner], cwd=PROJECT_ROOT)


def upload_wheels(wheels: list[Path], repository: str, assume_yes: bool) -> None:
    log("== 上传到 PyPI ==")
    for w in wheels:
        log(f"  - {w.name}")
    if not assume_yes:
        try:
            ans = input(f"确认把以上文件上传到 {repository}？[y/N] ").strip().lower()
        except EOFError:
            die("需要交互确认，但当前环境无法输入（可用 -y 跳过确认）")
        if ans not in ("y", "yes", "是"):
            log("已取消上传。")
            return
    twine = find_exe("twine")
    run([twine, "upload", "--repository", repository, *[str(w) for w in wheels]], cwd=PROJECT_ROOT)


# ------------------------------ 主流程 ------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(
        description="编译 Windows/Linux 两个平台的 wheel 并上传到 PyPI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-r", "--repository", default="pypi", help="目标仓库（pypi 或 testpypi）")
    parser.add_argument("--wsl-distro", default=DEFAULT_WSL_DISTRO, help="WSL 发行版名称")
    parser.add_argument("--only", choices=["both", "windows", "linux"], default="both", help="只编译某个平台")
    parser.add_argument("--no-upload", action="store_true", help="只编译，不上传")
    parser.add_argument("-y", "--yes", action="store_true", help="上传前不再二次确认")
    parser.add_argument("--skip-version-check", action="store_true", help="跳过与 PyPI 的版本冲突检查")
    args = parser.parse_args()

    name, version = read_project_meta()
    log(f"包名: {name}")
    log(f"本地版本: {version}")
    log(f"项目目录: {PROJECT_ROOT}")

    log("")
    if args.skip_version_check:
        log("[版本检查] 已通过 --skip-version-check 跳过。")
    else:
        check_version_conflict(name, version)

    log("")
    clean_dist()

    log("")
    if args.only in ("both", "windows"):
        build_windows(DIST_DIR)
    if args.only in ("both", "linux"):
        build_linux(DIST_DIR, args.wsl_distro)

    wheels = sorted(DIST_DIR.glob("*.whl"))
    if not wheels:
        die("没有在 dist/ 下找到任何生成的 wheel")

    log("")
    log("== 编译结果 ==")
    for w in wheels:
        log(f"  - {w.name}  ({w.stat().st_size / 1024 / 1024:.1f} MB)")

    if args.only == "both":
        if not any("win" in w.name for w in wheels):
            log("警告: 没有找到 Windows 平台的 wheel")
        if not any(("linux" in w.name or "manylinux" in w.name) for w in wheels):
            log("警告: 没有找到 Linux 平台的 wheel")

    log("")
    log("== 校验 wheel ==")
    run([find_exe("twine"), "check", *[str(w) for w in wheels]], cwd=PROJECT_ROOT)

    if args.no_upload:
        log("")
        log("已指定 --no-upload，跳过上传。")
        return

    log("")
    upload_wheels(wheels, args.repository, args.yes)
    log("")
    log("完成。")


if __name__ == "__main__":
    main()
