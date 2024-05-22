SMILES（Simplified Molecular Input Line Entry System）是一种用于表示化学分子结构的线性字符串表示法。生成SMILES表示的算法涉及遍历分子结构并将其转换为相应的字符串形式。以下是一个生成SMILES表示的简单算法概述：

1. **输入分子结构**：提供一个分子的图结构，通常是原子和键的邻接矩阵表示。

2. **选择起始原子**：选择一个原子作为起点。通常使用度数最低的原子作为起始点，如果有多个候选，则根据某种次序选择。

3. **深度优先搜索（DFS）**：
   - 从起始原子开始，执行深度优先搜索（DFS）。
   - 记录访问的路径，并跟踪已经访问的原子，以避免重复访问。

4. **生成部分字符串**：
   - 对于每个访问的原子，添加相应的原子符号（例如，C表示碳，O表示氧）。
   - 对于双键和三键，添加适当的符号（例如，=表示双键，#表示三键）。

5. **分支处理**：
   - 当遇到分支时，使用括号表示分支路径。
   - 每次遇到一个新的分支时，开启一个新的括号层，并在分支结束时关闭括号。

6. **环处理**：
   - 当遇到环结构时，使用数字标识环的起点和终点。
   - 在第一个遇到环的原子处添加一个数字，并在环闭合时添加相同的数字。

7. **输出SMILES字符串**：将所有部分字符串连接起来，生成最终的SMILES表示。

下面是一个Python伪代码示例，展示了如何实现简单的SMILES生成算法：

```python
class Atom:
    def __init__(self, symbol):
        self.symbol = symbol
        self.bonds = []

class Molecule:
    def __init__(self):
        self.atoms = []

def generate_smiles(molecule):
    visited = set()
    smiles = []
    stack = []

    def dfs(atom, parent=None):
        if atom in visited:
            return
        visited.add(atom)
        smiles.append(atom.symbol)
        for bond in atom.bonds:
            if bond != parent:
                smiles.append("-")  # Simple single bond representation
                dfs(bond, atom)
        if parent:
            smiles.append(")")

    for atom in molecule.atoms:
        if atom not in visited:
            smiles.append("(")
            dfs(atom)
    
    return "".join(smiles)

# Example usage
c1 = Atom("C")
c2 = Atom("C")
o = Atom("O")

c1.bonds.append(c2)
c2.bonds.append(c1)
c2.bonds.append(o)
o.bonds.append(c2)

molecule = Molecule()
molecule.atoms.extend([c1, c2, o])

smiles_string = generate_smiles(molecule)
print(smiles_string)  # Expected output: (CCO)
```

这个示例代码展示了如何通过深度优先搜索遍历分子结构并生成简单的SMILES字符串。当然，实际的SMILES生成算法会更加复杂，需要处理不同种类的键、原子和环结构，并且会考虑更复杂的化学规则。