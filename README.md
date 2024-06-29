# `pywfn` -- 基于python的波函数分析工具

详情请查看 https://www.xiaofei911.top/mkdocs/pywfn/


## 安装
1. 创建虚拟环境
```
python -m venv venv
```
2. 运行环境
```
.\venv\Scripts\activate
```
3. 安装依赖
```
pip install -r reqs.txt
```
4. 运行项目
```
python main.py
```


## 规范
- 传参的时候尽量不要使用自定义数据类型

- 当程序作为模块的时候，不要打印信息，尽量少使用print

代码的简化和重构也会花费大量时间呀

1.23,4.567

## 一些心得
- 一个模块中的所有文件中的类和函数都可以在__init__.py文件中引用
- 但是一个模块中的文件互相调用彼此的类的时候就不要使用__init__.py中的引用了
- 数据保存到py文件中，可以打包成单文件

## nuitka打包
> 自己的电脑打包有问题，可以用虚拟机来打包，使用Hyper-V

```
python -m nuitka --mingw64 --standalone --show-progress --output-dir=out main.py
python -m nuitka --mingw64 --standalone --show-progress --output-dir=out --windows-icon-from-ico=./ico.ico main.py
```

## 物理量缩写
- 轨道 obt
- 能量 eng
- 电子 ele elec
- 原子 atm atom
- 电荷 chg chrg

## 依赖的库
- numpy
- rich