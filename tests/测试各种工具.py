from pywfn.tools.log import Tool as LogTool

tool=LogTool()

opts=tool.split_opt(rf"d:\gfile\工具测试文件\lianxi_opt.out")
for opt in opts:
    atms,xyzs,eng=opt
    # print(atms)
    # print(xyzs)
    print(eng)
