import subprocess
import shutil

root='D:/code/pywfn'

result = subprocess.run(
    # ["cargo", "build", "--release"],
    ["cargo", "build"],
    cwd=root,  # 设置工作目录为 Rust 项目根目录
    capture_output=True,
    text=True
)

shutil.copy(f'{root}/target/debug/rlib.dll',f'{root}/pywfn/maths/rlib.pyd')
print(result)