use std::fs;
use std::path::PathBuf;

/// 构建脚本：从 Cargo.lock 中解析实际使用的 rswfn 版本号，
/// 通过 `cargo:rustc-env` 注入为编译期环境变量 `RSWFN_VERSION`，
/// 供 `src_rs` 中以 `env!("RSWFN_VERSION")` 读取。
fn main() {
    println!("cargo:rerun-if-changed=Cargo.lock");
    let version = rswfn_version().unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=RSWFN_VERSION={version}");
}

/// 在 Cargo.lock 中查找名为 `rswfn` 的包的版本号。
fn rswfn_version() -> Option<String> {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").ok()?;
    let lock = PathBuf::from(manifest_dir).join("Cargo.lock");
    let text = fs::read_to_string(lock).ok()?;

    // Cargo.lock 由若干 `[[package]]` 块组成，逐块解析 name/version
    for block in text.split("[[package]]") {
        let mut name = None;
        let mut version = None;
        for line in block.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("name = ") {
                name = Some(rest.trim().trim_matches('"').to_string());
            } else if let Some(rest) = line.strip_prefix("version = ") {
                version = Some(rest.trim().trim_matches('"').to_string());
            }
        }
        if name.as_deref() == Some("rswfn") {
            return version;
        }
    }
    None
}
