"""
将需要的项目代码打包
"""
from pathlib import Path
import zipfile
import os

class Tool:
    def __init__(self) -> None:
        self.forders=[
            'pywfn',
        ]
        self.files=[
            'main.py',
            'main.bat',
            'reqs.txt'
        ]
        # 定义需要排除的文件类型
        self.rms=[
            '__pycache__',
            '.pyc',
            '.ui',
        ]

    def isRemove(self,name):
        for rm in self.rms:
            if rm in name:
                return True
        return False


    def zip(self,name):
        for forder in self.forders:
            self.files+=self.get_files(forder)
        for i,file in enumerate(self.files):
            print(i,file)
        zipFile=zipfile.ZipFile(f'{name}.zip','w')
        for each in self.files:
            zipFile.write(each)
        
    def get_files(self,path):
        allFiles = []
        for root, dirs, files in os.walk(path):
            for file in files:
                name=os.path.join(root, file)
                if self.isRemove(name):
                    continue
                allFiles.append(name)
        return allFiles
    
tool=Tool()
tool.zip('pywfn')