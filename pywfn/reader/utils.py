import numpy as np

def toCart(atms:list[int],shls:list[int],syms:list[str],CM:np.ndarray): # 将数据转为笛卡尔类型的
    CMlist=[]
    atmList=[]
    shlList=[]
    symList=[]
    i=0
    from pywfn.data import bastrans as bt
    while i<CM.shape[0]:
        sym=syms[i]
        match sym:
            case 'D 0':
                CMlist.append(bt.DMat@CM[i:i+5])
                atmList+=[atms[i]]*6
                shlList+=[shls[i]]*6
                symList+=bt.Dsyms
                i+=5
            case 'F 0':
                CMlist.append(bt.FMat@CM[i:i+7])
                atmList+=[atms[i]]*10
                shlList+=[shls[i]]*10
                symList+=bt.Fsyms
                i+=7
            case 'G 0':
                CMlist.append(bt.GMat@CM[i:i+9])
                atmList+=[atms[i]]*15
                shlList+=[shls[i]]*15
                symList+=bt.Gsyms
                i+=9
            case 'H 0':
                CMlist.append(bt.HMat@CM[i:i+11])
                atmList+=[atms[i]]*21
                shlList+=[shls[i]]*21
                symList+=bt.Hsyms
                i+=11
            case _:
                CMlist.append(CM[i][None,:])
                atmList.append(atms[i])
                shlList.append(shls[i])
                symList.append(syms[i])
                i+=1

    CM=np.concatenate(CMlist,axis=0)
    return atmList,shlList,symList,CM