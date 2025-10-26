将当前的模块添加到父模块之中

```rust
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "reads")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(函数名, &m)?)?; // 添加这个rs文件中的函数到子模块中
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
```
