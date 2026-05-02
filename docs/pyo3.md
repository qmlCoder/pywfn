## 矩阵的项目转换

rust中使用的矩阵为ndarray，python中使用的矩阵为numpy

rust接收python矩阵

```rust
def fn set_cmat(arr:cmat: Bound<'_, PyArray2<f64>>){
    let cmat=cmat.into_array();
}
```

rust发送python矩阵

```rust
def fn get_cmat(py:Python)->Py<PyArray2<f64>>{
    cmat.into_pyarray(py).unbind()
}
```
