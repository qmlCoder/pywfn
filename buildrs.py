import subprocess
import shutil

root='D:/code/pywfn'

release=False

if release:
    result = subprocess.run(
        ["cargo", "build", "--release"],
        cwd=root,  # 设置工作目录为 Rust 项目根目录
        capture_output=True,
        text=True
    )
    shutil.copy(f'{root}/target/release/core.dll',f'{root}/pywfn/core.pyd')
else:
    result = subprocess.run(
        ["cargo", "build"],
        cwd=root,  # 设置工作目录为 Rust 项目根目录
        capture_output=True,
        text=True
    )
    shutil.copy(f'{root}/target/debug/core.dll',f'{root}/pywfn/core.pyd')
print(result)