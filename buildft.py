import subprocess
import shutil

root='D:/code/pywfn'

result = subprocess.run(
    # ["cargo", "build", "--release"],
    ["gfortran", "-shared", "*.f90", "-o", "flib.dll"],
    cwd=f'{root}/libft',  # 设置工作目录为 Rust 项目根目录
    capture_output=True,
    text=True
)

# shutil.copy(f'{root}/libft/flib.dll',f'{root}/pywfn/maths/flib.dll')
# print(result)