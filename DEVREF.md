# 开发参考

核心代码逻辑在rswfn中完成，本项目现在是套壳

套了两层壳：
- 第一层壳是rust代码生成的python二进制模块：core.pyd
- 第二层壳是python代码对core.pyd的封装，方便有类型提示

能使用numpy数组的都要使用numpy数组

直接使用pyi文件，不写自己的代码逻辑了