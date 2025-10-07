# Markdown渲染测试

这是一个测试文件，用于验证markdown渲染功能。

## 代码块测试

### Python代码
```python
def hello_world():
    print("Hello, World!")
    return "success"

# 这是一个注释
x = [1, 2, 3, 4, 5]
for i in x:
    if i % 2 == 0:
        print(f"{i} 是偶数")
```

### JavaScript代码
```javascript
function fibonacci(n) {
    if (n <= 1) return n;
    return fibonacci(n - 1) + fibonacci(n - 2);
}

// 数组方法示例
const numbers = [1, 2, 3, 4, 5];
const doubled = numbers.map(x => x * 2);
console.log(doubled); // [2, 4, 6, 8, 10]
```

### Bash命令
```bash
# 列出文件
ls -la

# 创建目录
mkdir test_dir
cd test_dir
```

## HTML支持测试

<div style="background: #f0f0f0; padding: 10px; border-radius: 5px;">
  <h3 style="color: #333;">自定义HTML内容</h3>
  <p>这是一个HTML段落，测试HTML标签的支持。</p>
</div>

## 表格测试

| 姓名 | 年龄 | 职业 | 城市 |
|------|------|------|------|
| 张三 | 25 | 工程师 | 北京 |
| 李四 | 30 | 设计师 | 上海 |
| 王五 | 28 | 产品经理 | 深圳 |

## 引用测试

> 这是一个引用示例。
>
> 支持多行引用内容。
>
> > 嵌套引用也是支持的。

## 列表测试

### 无序列表
- 项目1
- 项目2
  - 子项目2.1
  - 子项目2.2
- 项目3

### 有序列表
1. 第一步
2. 第二步
3. 第三步

## 链接和图片

[这是一个链接](https://example.com)

![图片示例](https://via.placeholder.com/300x200)

## 行内代码

这里有一些行内代码：`const x = 10;` 和 `print("hello")`。

## 分割线

---

## 文本样式

**粗体文本**、*斜体文本*、***粗斜体文本***、~~删除线~~

## 数学公式（如果支持）

$E = mc^2$

$$
\sum_{i=1}^{n} i = \frac{n(n+1)}{2}
$$