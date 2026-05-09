<h1 align = "center">矩阵理论与应用II习题</h1>
<p align="right">SY2502110-朱宛瑜<br>Date: 2026.5.6</p>

# 1.自测题

1. 设 $A, B \in \mathbb{C}^{n \times n}$，且 $A$ 的特征值全为实数，证明：矩阵方程 $X + AXA + A^2XA^2 = B$ 有唯一解。

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **定理1.2.5（直积/列拉直）**（见《矩阵理论与应用II复习》 §1.1.2 方法三）与 **推论1.2.6**（直积特征值性质）。将矩阵方程拉直为线性方程组 $M \cdot \text{vec}(X) = \text{vec}(B)$，证明系数矩阵 $M$ 可逆即可。
   >
   > **【证明】**
   >
   > 对原方程两边列拉直，利用定理1.2.5：$vec(AXB) = (B^T \otimes A)vec(X)$，有：
   > $$
   > vec(X) + vec(AXA) + vec(A^2XA^2) = vec(B)\\
   > vec(X) + (A^T \otimes A)vec(X) + ((A^2)^T \otimes A^2)vec(X) = vec(B)
   > $$
   > 即：
   > $$
   > (I_{n^2} + A^T \otimes A + (A^T)^2 \otimes A^2) \cdot vec(X) = vec(B)
   > $$
   > 记 $M = I + A^T \otimes A + (A^T)^2 \otimes A^2 \in \mathbb{C}^{n^2 \times n^2}$。
   >
   > 方程有唯一解当且仅当 $M$ 可逆（即 $\det(M) \neq 0$），等价于 $M$ 无零特征值。
   >
   > 由推论1.2.6，$A^T \otimes A$ 的特征值为 $\lambda_i \lambda_j$（$i,j=1,\ldots,n$），其中 $\lambda_i$ 为 $A$ 的特征值。同理$(A^T)^2 \otimes A^2$ 的特征值为 $\lambda_i^2 \lambda_j^2$。
   >
   > 故 $M$ 的特征值为 $\mu_{ij} = 1 + \lambda_i \lambda_j + \lambda_i^2 \lambda_j^2$（由矩阵特征值的性质所得，矩阵多项式的特征值为将特征值代入多项式表达式所得到的值）。
   >
   > 因为 ==$A$ 的特征值全为实数==，令 $t = \lambda_i \lambda_j \in \mathbb{R}$，则 $\mu_{ij} = 1 + t + t^2$。二次函数 $f(t) = t^2 + t + 1$ 的判别式 $\Delta = 1 - 4 = -3 < 0$，且 $f(0) = 1 > 0$，故对任意实数 $t$，恒有 $f(t) > 0$。因此 $M$ 的所有特征值均不为零，$M$ 可逆，原矩阵方程有唯一解。
   >

2. 证明：$A \in \mathbb{R}_r^{m \times n}$，$\sigma_i$ ($1 \leqslant i \leqslant r$) 为奇异值，则 $\text{tr}(A^H A) = \sum_{i=1}^r \sigma_i^2$。

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **定义1.1.2（奇异值）** 与 **引理1.1.3**（见《矩阵理论与应用II复习》 §1.1.1）：矩阵 $A^HA$ 的非零特征值个数为 $r$，奇异值 $\sigma_i = \sqrt{\lambda_i}$，其中 $\lambda_i$ 是 $A^HA$ 的特征值。矩阵的迹等于其特征值之和。
   >
   > **【证明】**
   >
   > 由定义1.1.2，设 $A^HA \in \mathbb{R}^{n \times n}$ 的特征值为 $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_r > 0$，$\lambda_{r+1} = \cdots = \lambda_n = 0$，则 $A$ 的正奇异值为 $\sigma_i = \sqrt{\lambda_i}$，$i=1,\ldots,r$。
   >
   > 矩阵的迹等于其特征值之和，故：
   > $$
   > \text{tr}(A^HA) = \sum_{i=1}^n \lambda_i = \sum_{i=1}^r \lambda_i + \sum_{i=r+1}^n 0 = \sum_{i=1}^r \lambda_i
   > $$
   > 代入 $\lambda_i = \sigma_i^2$，即得：
   > $$
   > \text{tr}(A^HA) = \sum_{i=1}^r \sigma_i^2
   > $$

3. 求下列矩阵的奇异值分解：
   $$
   (1) \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix},
   \quad (2) \begin{pmatrix} 1 & -1 \\ 3 & -3 \\ -3 & 3 \end{pmatrix},
   \quad (3) \begin{pmatrix} 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix},
   \quad (4) \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}.
   $$

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **定理1.1.1（奇异值分解）** 与 **奇异值分解步骤2**（见《矩阵理论与应用II复习》 §1.1.1）：先计算 $A^HA$ 的特征值与标准正交特征向量构成 $V$，再由 $u_i = \frac{1}{\sigma_i} A v_i$ 求 $U$ 的前 $r$ 列并扩充为标准正交基。优先计算阶数较小的矩阵。
   >
   > ---
   >
   > **（1）** $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix}$（$3 \times 2$）
   >
   > **步骤一：** 计算 $A^HA$（$2 \times 2$，阶数更小）：
   > $$
   > A^HA = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}
   > $$
   > **步骤二：** 求 $A^HA$ 的特征值：
   > $$
   > |\lambda I - A^HA| = \begin{vmatrix} \lambda-2 & -1 \\ -1 & \lambda-2 \end{vmatrix} = (\lambda-2)^2 - 1 = (\lambda-1)(\lambda-3) = 0
   > $$
   > 得 $\lambda_1 = 3, \lambda_2 = 1$。奇异值 $\sigma_1 = \sqrt{3}, \sigma_2 = 1$，$r = 2$。
   >
   > **步骤三：** 求标准正交特征向量构成 $V$：
   >
   > - $\lambda_1 = 3$：$(3I - A^HA)x = 0 \Rightarrow \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}x = 0 \Rightarrow x_1 = x_2$，取 $v_1 = \frac{1}{\sqrt{2}}(1, 1)^T$；
   > - $\lambda_2 = 1$：$(I - A^HA)x = 0 \Rightarrow \begin{pmatrix} -1 & -1 \\ -1 & -1 \end{pmatrix}x = 0 \Rightarrow x_1 = -x_2$，取 $v_2 = \frac{1}{\sqrt{2}}(1, -1)^T$，得到
   >
   > $$
   > V = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}
   > $$
   >
   > **步骤四：** 求 $U$ 的前 $r$ 列：
   > $$
   > u_1 = \frac{1}{\sigma_1} A v_1 = \frac{1}{\sqrt{3}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix} \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{1}{\sqrt{6}}\begin{pmatrix} 1 \\ 1 \\ 2 \end{pmatrix}
   > $$
   >
   > $$
   > u_2 = \frac{1}{\sigma_2} A v_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 1 \end{pmatrix} \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ -1 \\ 0 \end{pmatrix}
   > $$
   >
   > **步骤五：** 扩充 $u_3$ 为标准正交基。设 $u_3 = (x, y, z)^T$，满足 $u_1^T u_3 = 0, u_2^T u_3 = 0$：
   > $$
   > \begin{cases} x + y + 2z = 0 \\ x - y = 0 \end{cases} \Rightarrow x = y,\; z = -x
   > $$
   > 取 $u_3 = \frac{1}{\sqrt{3}}(1, 1, -1)^T$，我们有
   > $$
   > U = \begin{pmatrix} \frac{1}{\sqrt{6}} & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{3}} \\[4pt] \frac{1}{\sqrt{6}} & -\frac{1}{\sqrt{2}} & \frac{1}{\sqrt{3}} \\[4pt] \frac{2}{\sqrt{6}} & 0 & -\frac{1}{\sqrt{3}} \end{pmatrix}
   > $$
   > **奇异值分解：** $A = U \begin{pmatrix} \sqrt{3} & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} V^H$
   >
   > ---
   >
   > **（2）** $A = \begin{pmatrix} 1 & -1 \\ 3 & -3 \\ -3 & 3 \end{pmatrix}$（$3 \times 2$）
   >
   > 观察：$A = \begin{pmatrix} 1 \\ 3 \\ -3 \end{pmatrix} \begin{pmatrix} 1 & -1 \end{pmatrix}$，为秩1矩阵。
   >
   > **步骤一：** $A^HA = \begin{pmatrix} 1 & 3 & -3 \\ -1 & -3 & 3 \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 3 & -3 \\ -3 & 3 \end{pmatrix} = \begin{pmatrix} 19 & -19 \\ -19 & 19 \end{pmatrix}$
   >
   > **步骤二：** $|\lambda I - A^HA| = \begin{vmatrix} \lambda-19 & 19 \\ 19 & \lambda-19 \end{vmatrix} = (\lambda-19)^2 - 361 = \lambda(\lambda-38) = 0$
   >
   > $\lambda_1 = 38, \lambda_2 = 0$。$\sigma_1 = \sqrt{38}$，$r = 1$。
   >
   > **步骤三：** 特征向量（已正交）：
   > - $\lambda_1 = 38$：$\begin{pmatrix} 19 & 19 \\ 19 & 19 \end{pmatrix}x = 0 \Rightarrow x_1 = -x_2$，$v_1 = \frac{1}{\sqrt{2}}(1, -1)^T$；
   > - $\lambda_2 = 0$：$v_2 = \frac{1}{\sqrt{2}}(1, 1)^T$，得到
   >
   > $$
   > V = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}
   > $$
   >
   > **步骤四：** $u_1 = \frac{1}{\sigma_1} A v_1 = \frac{1}{\sqrt{38}} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix} 2 \\ 6 \\ -6 \end{pmatrix} = \frac{1}{\sqrt{19}}\begin{pmatrix} 1 \\ 3 \\ -3 \end{pmatrix}$
   >
   > **步骤五：** 扩充 $u_2, u_3$。由 $u_2 \perp u_1$：取 $u_2 = \frac{1}{\sqrt{10}}(-3, 1, 0)^T$。由 $u_3 \perp u_1, u_2$（叉积）：$u_3 = \frac{1}{\sqrt{190}}(3, 9, 10)^T$，将三个 $u$ 按照对应顺序拼起来，得到
   > $$
   > U = \begin{pmatrix} \frac{1}{\sqrt{19}} & -\frac{3}{\sqrt{10}} & \frac{3}{\sqrt{190}} \\[4pt] \frac{3}{\sqrt{19}} & \frac{1}{\sqrt{10}} & \frac{9}{\sqrt{190}} \\[4pt] -\frac{3}{\sqrt{19}} & 0 & \frac{10}{\sqrt{190}} \end{pmatrix}
   > $$
   > **奇异值分解：** $A = U \begin{pmatrix} \sqrt{38} & 0 \\ 0 & 0 \\ 0 & 0 \end{pmatrix} V^H$
   >
   > ---
   >
   > **（3）** $A = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$（$2 \times 3$，优先计算 $AA^H$）
   >
   > **步骤一：** $AA^H = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 1 & 0 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$
   >
   > **步骤二：** 特征值 $2, 1$。$\sigma_1 = \sqrt{2}, \sigma_2 = 1$，$r = 2$。$U = I_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$。
   >
   > **步骤三：** 由 $V_1 = A^H U_1 \Sigma_r^{-1}$：
   > $$
   > v_1 = \frac{A^H u_1}{\sigma_1} = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \\ 0 \end{pmatrix}, \quad v_2 = \frac{A^H u_2}{\sigma_2} = \begin{pmatrix} 0 \\ 0 \\ 1 \end{pmatrix}
   > $$
   >
   > **步骤四：** 扩充 $v_3 \perp v_1, v_2$，得到$v_3 = \frac{1}{\sqrt{2}}(1, -1, 0)^T$：
   > $$
   > V = \begin{pmatrix} \frac{1}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} \\[4pt] \frac{1}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}} \\[4pt] 0 & 1 & 0 \end{pmatrix}
   > $$
   > **奇异值分解：** $A = \begin{pmatrix} \sqrt{2} & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} V^H$（$U = I_2$）
   >
   > ---
   >
   > **（4）** $A = \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$（$3 \times 1$）
   >
   > **步骤一：** $A^HA = \begin{pmatrix} 1 & 1 & 1 \end{pmatrix} \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} = 3$。$\lambda_1 = 3, \sigma_1 = \sqrt{3}, r = 1$。
   >
   > **步骤二：** $v_1 = 1$（$V$ 为 $1 \times 1$ 单位阵）。
   >
   > **步骤三：** $u_1 = \frac{1}{\sigma_1} A v_1 = \frac{1}{\sqrt{3}}\begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix}$。
   >
   > **步骤四：** 扩充 $u_2 = \frac{1}{\sqrt{2}}(1, -1, 0)^T$，$u_3 = u_1 \times u_2 = \frac{1}{\sqrt{6}}(1, 1, -2)^T$，得到：
   > $$
   > U = \begin{pmatrix} \frac{1}{\sqrt{3}} & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{6}} \\[4pt] \frac{1}{\sqrt{3}} & -\frac{1}{\sqrt{2}} & \frac{1}{\sqrt{6}} \\[4pt] \frac{1}{\sqrt{3}} & 0 & -\frac{2}{\sqrt{6}} \end{pmatrix}
   > $$
   > **奇异值分解：** $A = U \begin{pmatrix} \sqrt{3} \\ 0 \\ 0 \end{pmatrix} \cdot 1$

4. 设 $A$ 是正规矩阵，证明：
   （1） $A^+ A = AA^+$

   （2） $(A^+)^k = (A^k)^+$ （$k$为正整数）

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **命题1.1.1（正规矩阵的奇异值）** 以及正规矩阵可酉对角化的性质（见《矩阵理论与应用II复习》 §1.1.1 与《矩阵论复习》第17条）。正规矩阵 $A$ 满足 $A^HA = AA^H$，必可酉对角化：$A = U \Lambda U^H$，其中 $U$ 为酉矩阵，$\Lambda = \text{diag}(\lambda_1, \ldots, \lambda_n)$。再由 **定理1.6.1** 中加号逆的表达式，$A^+ = U \Lambda^+ U^H$，其中 $\lambda_i^+ = \begin{cases} \lambda_i^{-1} & \lambda_i \neq 0 \\ 0 & \lambda_i = 0 \end{cases}$。
   >
   > **（1）证明 $A^+A = AA^+$：**
   >
   > 由酉对角化 $A = U \Lambda U^H$，$A^+ = U \Lambda^+ U^H$，有：
   > $$
   > A^+A = U \Lambda^+ U^H \cdot U \Lambda U^H = U \Lambda^+ \Lambda U^H\\
   > AA^+ = U \Lambda U^H \cdot U \Lambda^+ U^H = U \Lambda \Lambda^+ U^H
   > $$
   > 对角矩阵 $\Lambda$ 与 $\Lambda^+$ 可交换（$\Lambda^+\Lambda = \Lambda\Lambda^+$ 对角元为 $0$ 或 $1$），故 $A^+A = AA^+$。
   >
   > **（2）证明 $(A^+)^k = (A^k)^+$：**
   > $$
   > A^k = (U \Lambda U^H)^k = U \Lambda^k U^H
   > $$
   > 由加号逆定义：
   > $$
   > (A^k)^+ = U (\Lambda^k)^+ U^H
   > $$
   > 而 $(\Lambda^k)^+ = \text{diag}((\lambda_1^k)^+, \ldots, (\lambda_n^k)^+)$，且 $(\lambda_i^k)^+ = (\lambda_i^+)^k$（验证：若 $\lambda_i \neq 0$，则 $(\lambda_i^k)^{-1} = (\lambda_i^{-1})^k$；若 $\lambda_i = 0$，两者均为 $0$）。故：
   > $$
   > (A^k)^+ = U (\Lambda^+)^k U^H = (U \Lambda^+ U^H)^k = (A^+)^k
   > $$

5. 验证方程是相容的，并用 $A^{(1)}$ 表示通解。
   $$
   \begin{cases}
   x_1 + 2x_2 - x_3 = 1 \\
   -x_2 + 2x_3 = 2
   \end{cases}
   $$

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **推论1.2.3（Penrose定理）** 与 **命题1.3.2（相容条件）**（见《矩阵理论与应用II复习》 §1.1.2 方法二 及 §1.1.3）。对相容方程组 $Ax = b$，通解为 $x = A^{(1)}b + (I - A^{(1)}A)y$（$\forall y \in \mathbb{C}^n$）。
   >
   > **步骤一：验证相容性。**
   >
   > $A = \begin{pmatrix} 1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix},b = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$，我们有$\text{rank}(A) = 2$（两行线性无关）。增广矩阵：
   > $$
   > (A, b) = \begin{pmatrix} 1 & 2 & -1 & 1 \\ 0 & -1 & 2 & 2 \end{pmatrix},\; \text{rank}(A, b) = 2
   > $$
   > $\text{rank}(A) = \text{rank}(A, b) = 2$，故方程相容（命题1.3.2）。
   >
   > **步骤二：求 $A^{(1)}$。** $A$ 为行满秩矩阵（$2 \times 3$，秩 $2$），使用右逆（《矩阵理论与应用II复习》 §1.1.3 第1点）：
   > $$
   > A^{(1)} = A^T(AA^T)^{-1}
   > $$
   >
   > $$
   > AA^T = \begin{pmatrix} 1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 2 & -1 \\ -1 & 2 \end{pmatrix} = \begin{pmatrix} 6 & -4 \\ -4 & 5 \end{pmatrix}
   > $$
   > $$
   > (AA^T)^{-1} = \frac{1}{30-16}\begin{pmatrix} 5 & 4 \\ 4 & 6 \end{pmatrix} = \frac{1}{14}\begin{pmatrix} 5 & 4 \\ 4 & 6 \end{pmatrix}
   > $$
   > $$
   > A^{(1)} = \frac{1}{14}\begin{pmatrix} 1 & 0 \\ 2 & -1 \\ -1 & 2 \end{pmatrix} \begin{pmatrix} 5 & 4 \\ 4 & 6 \end{pmatrix} = \frac{1}{14}\begin{pmatrix} 5 & 4 \\ 6 & 2 \\ 3 & 8 \end{pmatrix}
   > $$
   >
   > 验证：$AA^{(1)} = I_2$，确为减号逆。
   >
   > **步骤三：写出通解。** 由推论1.2.3：
   > $$
   > x = A^{(1)}b + (I - A^{(1)}A)y, \quad \forall y = (y_1, y_2, y_3)^T \in \mathbb{C}^3
   > $$
   >
   > $$
   > A^{(1)}b = \frac{1}{14}\begin{pmatrix} 5 & 4 \\ 6 & 2 \\ 3 & 8 \end{pmatrix}\begin{pmatrix} 1 \\ 2 \end{pmatrix} = \frac{1}{14}\begin{pmatrix} 13 \\ 10 \\ 19 \end{pmatrix}
   > $$
   >
   > $$
   > A^{(1)}A = \frac{1}{14}\begin{pmatrix} 5 & 4 \\ 6 & 2 \\ 3 & 8 \end{pmatrix}\begin{pmatrix} 1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix} = \frac{1}{14}\begin{pmatrix} 5 & 6 & 3 \\ 6 & 10 & -2 \\ 3 & -2 & 13 \end{pmatrix}
   > $$
   >
   > $$
   > I - A^{(1)}A = \frac{1}{14}\begin{pmatrix} 9 & -6 & -3 \\ -6 & 4 & 2 \\ -3 & 2 & 1 \end{pmatrix}
   > $$
   >
   > 注意到 $(I - A^{(1)}A)y = \frac{3y_1 - 2y_2 - y_3}{14}\begin{pmatrix} 3 \\ -2 \\ -1 \end{pmatrix}$。令 $t = \frac{3y_1 - 2y_2 - y_3}{14}$，通解亦可简洁表为：
   > $$
   > x = \frac{1}{14}\begin{pmatrix} 13 \\ 10 \\ 19 \end{pmatrix} + t\begin{pmatrix} 3 \\ -2 \\ -1 \end{pmatrix}, \quad \forall t \in \mathbb{C}
   > $$

6. 求上题通解中的极小范数解。

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **定理1.4.4** 与 **推论1.4.1**（见《矩阵理论与应用II复习》 §1.1.4）：$G \in A\{1,4\}$ 的充要条件是 $x = Gb$ 为相容方程组的极小范数解，且极小范数解唯一。对于行满秩矩阵，$A^H(AA^H)^{-1} \in A\{1,4\}$。等价地，由 **定理1.6.5**（§1.1.6），$x = A^+b$ 是相容方程的唯一极小范数解。
   >
   > 上题中 $A$ 行满秩，$A^+ = A^T(AA^T)^{-1}$（推论1.6.2），恰与第5题中构造的 $A^{(1)}$ 相同。故极小范数解为：
   > $$
   > x_{\min} = A^+b = \frac{1}{14}\begin{pmatrix} 13 \\ 10 \\ 19 \end{pmatrix}
   > $$
   >
   > 对应通解中 $t = 0$ 的特解。由推论1.4.1，极小范数解唯一。

7. 求解矩阵方程 $AX + XB = C$，其中

   （1） $A = \begin{pmatrix} 1 & -1 \\ 0 & 2 \end{pmatrix}, B = \begin{pmatrix} -3 & 4 \\ 1 & 0 \end{pmatrix}, C = \begin{pmatrix} 1 & 3 \\ -2 & 2 \end{pmatrix}$；

   （2） $A = \begin{pmatrix} 1 & -1 \\ 0 & 2 \end{pmatrix}, B = \begin{pmatrix} -3 & 4 \\ 0 & -1 \end{pmatrix}, C = \begin{pmatrix} 0 & 5 \\ 2 & -9 \end{pmatrix}$。

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **方法三：直积和矩阵拉直**（见《矩阵理论与应用II复习》 §1.1.2 方法三）。由 **定理1.2.5**：$vec(AX) = (I \otimes A)vec(X)$，$vec(XB) = (B^T \otimes I)vec(X)$。方程 $AX + XB = C$ 等价于：
   > $$
   > (I \otimes A + B^T \otimes I) \cdot vec(X) = vec(C)
   > $$
   > 记 $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$，$vec(X) = (x_{11}, x_{21}, x_{12}, x_{22})^T$。
   >
   > ---
   >
   > **（1）** $A = \begin{pmatrix} 1 & -1 \\ 0 & 2 \end{pmatrix}, B = \begin{pmatrix} -3 & 4 \\ 1 & 0 \end{pmatrix}$
   >
   > $$
   > I \otimes A = \begin{pmatrix} A & 0 \\ 0 & A \end{pmatrix} = \begin{pmatrix} 1 & -1 & 0 & 0 \\ 0 & 2 & 0 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 2 \end{pmatrix}
   > $$
   >
   > $$
   > B^T = \begin{pmatrix} -3 & 1 \\ 4 & 0 \end{pmatrix},\; B^T \otimes I = \begin{pmatrix} -3I & I \\ 4I & 0 \end{pmatrix} = \begin{pmatrix} -3 & 0 & 1 & 0 \\ 0 & -3 & 0 & 1 \\ 4 & 0 & 0 & 0 \\ 0 & 4 & 0 & 0 \end{pmatrix}
   > $$
   >
   > 系数矩阵 $M = I \otimes A + B^T \otimes I$：
   > $$
   > M = \begin{pmatrix} -2 & -1 & 1 & 0 \\ 0 & -1 & 0 & 1 \\ 4 & 0 & 1 & -1 \\ 0 & 4 & 0 & 2 \end{pmatrix},\quad vec(C) = \begin{pmatrix} 1 \\ -2 \\ 3 \\ 2 \end{pmatrix}
   > $$
   >
   > 解 $M \cdot vec(X) = vec(C)$：
   > - 第2行：$-x_{21} + x_{22} = -2 \Rightarrow x_{22} = x_{21} - 2$
   > - 第4行：$4x_{21} + 2x_{22} = 2 \Rightarrow 4x_{21} + 2(x_{21}-2) = 2 \Rightarrow 6x_{21} = 6 \Rightarrow x_{21} = 1$，故 $x_{22} = -1$
   > - 第1行：$-2x_{11} - 1 + x_{12} = 1 \Rightarrow -2x_{11} + x_{12} = 2$
   > - 第3行：$4x_{11} + x_{12} - (-1) = 3 \Rightarrow 4x_{11} + x_{12} = 2$
   >
   > 两式相减：$(4x_{11} + x_{12}) - (-2x_{11} + x_{12}) = 2 - 2 \Rightarrow 6x_{11} = 0 \Rightarrow x_{11} = 0$，故 $x_{12} = 2$。
   >
   > $$
   > X = \begin{pmatrix} 0 & 2 \\ 1 & -1 \end{pmatrix}
   > $$
   >
   > 验证：$AX + XB = \begin{pmatrix} -1 & 3 \\ 2 & -2 \end{pmatrix} + \begin{pmatrix} 2 & 0 \\ -4 & 4 \end{pmatrix} = \begin{pmatrix} 1 & 3 \\ -2 & 2 \end{pmatrix} = C$。✓
   >
   > ---
   >
   > **（2）** $A = \begin{pmatrix} 1 & -1 \\ 0 & 2 \end{pmatrix}, B = \begin{pmatrix} -3 & 4 \\ 0 & -1 \end{pmatrix}$
   >$$
   > B^T = \begin{pmatrix} -3 & 0 \\ 4 & -1 \end{pmatrix},\; B^T \otimes I = \begin{pmatrix} -3I & 0 \\ 4I & -I \end{pmatrix} = \begin{pmatrix} -3 & 0 & 0 & 0 \\ 0 & -3 & 0 & 0 \\ 4 & 0 & -1 & 0 \\ 0 & 4 & 0 & -1 \end{pmatrix}
   > $$
   > 
   >$$
   > M = \begin{pmatrix} -2 & -1 & 0 & 0 \\ 0 & -1 & 0 & 0 \\ 4 & 0 & 0 & -1 \\ 0 & 4 & 0 & 1 \end{pmatrix},\quad vec(C) = \begin{pmatrix} 0 \\ 2 \\ 5 \\ -9 \end{pmatrix}
   > $$
   > 
   >解方程组：
   > - 第2行：$-x_{21} = 2 \Rightarrow x_{21} = -2$
   > - 第1行：$-2x_{11} - x_{21} = 0 \Rightarrow -2x_{11} + 2 = 0 \Rightarrow x_{11} = 1$
   > - 第3行：$4x_{11} - x_{22} = 5 \Rightarrow 4(1) - x_{22} = 5 \Rightarrow x_{22} = -1$
   > - 第4行：$4x_{21} + x_{22} = -9 \Rightarrow 4(-2) + (-1) = -9$ ✓（恒成立）
   > 
   >$x_{12}$ 未出现在方程中，取值任意。故：
   > $$
   > X = \begin{pmatrix} 1 & x_{12} \\ -2 & -1 \end{pmatrix}, \quad \forall x_{12} \in \mathbb{C}
   > $$
   > 
   >验证：$AX + XB = \begin{pmatrix} 3 & x_{12}+1 \\ -4 & -2 \end{pmatrix} + \begin{pmatrix} -3 & 4-x_{12} \\ 6 & -7 \end{pmatrix} = \begin{pmatrix} 0 & 5 \\ 2 & -9 \end{pmatrix} = C$。✓
   
8. 设 $A = \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 2 & 2 & 2 \\ -1 & 4 & 5 & 3 \end{pmatrix}, \mathbf{b} = \begin{pmatrix} 4 \\ -2 \\ -2 \end{pmatrix}$，
   （1）用广义逆矩阵方法判断方程组 $A\mathbf{x} = \mathbf{b}$ 是否相容；
   （2）求 $A\mathbf{x} = \mathbf{b}$ 的极小范数解或极小最小二乘解。

   > [!IMPORTANT]
   >
   > **【原理依据】**
   > - 相容性判断：**命题1.3.2** 与 Penrose定理（《矩阵理论与应用II复习》 §1.1.3）：$AA^{(1)}b = b$ 当且仅当方程组相容。或直接用 $AA^+b = b$（推论1.6.3）。
   > - 求解 $A^+$：**定理1.6.4**（满秩分解求加号逆）（《矩阵理论与应用II复习》 §1.1.6）：$A^+ = C^H(CC^H)^{-1}(B^HB)^{-1}B^H$。
   > - 最佳逼近解：**定理1.6.5**（《矩阵理论与应用II复习》 §1.1.6）：对不相容方程组，$x = A^+b$ 是唯一的极小最小二乘解（唯一最佳逼近解）；对相容方程组则是唯一的极小范数解。
   >
   > ---
   >
   > **步骤一：满秩分解。**
   >
   > 对 $A$ 做初等行变换化为行阶梯形：
   > $$
   > \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 2 & 2 & 2 \\ -1 & 4 & 5 & 3 \end{pmatrix} \xrightarrow{R_3 + R_1} \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 2 & 2 & 2 \\ 0 & 4 & 4 & 4 \end{pmatrix} \xrightarrow{R_3 - 2R_2} \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 2 & 2 & 2 \\ 0 & 0 & 0 & 0 \end{pmatrix}
   > $$
   >
   > $\text{rank}(A) = 2$。取前两列为列满秩矩阵 $B$：
   > $$
   > B = \begin{pmatrix} 1 & 0 \\ 0 & 2 \\ -1 & 4 \end{pmatrix}
   > $$
   >
   > 将 $A$ 的每一列用 $B$ 的列线性表出：
   > - $a_1 = 1 \cdot b_1 + 0 \cdot b_2 \Rightarrow c_1 = (1, 0)^T$
   > - $a_2 = 0 \cdot b_1 + 1 \cdot b_2 \Rightarrow c_2 = (0, 1)^T$
   > - $a_3 = -1 \cdot b_1 + 1 \cdot b_2 \Rightarrow c_3 = (-1, 1)^T$
   > - $a_4 = 1 \cdot b_1 + 1 \cdot b_2 \Rightarrow c_4 = (1, 1)^T$
   >
   > $$
   > C = \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 1 & 1 & 1 \end{pmatrix},\quad A = BC
   > $$
   >
   > **步骤二：用定理1.6.4 求 $A^+$。**
   >
   > $$
   > CC^H = \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 1 & 1 & 1 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ -1 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 3 & 0 \\ 0 & 3 \end{pmatrix} = 3I,\quad (CC^H)^{-1} = \frac{1}{3}I
   > $$
   >
   > $$
   > B^HB = \begin{pmatrix} 1 & 0 & -1 \\ 0 & 2 & 4 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 2 \\ -1 & 4 \end{pmatrix} = \begin{pmatrix} 2 & -4 \\ -4 & 20 \end{pmatrix}
   > $$
   > $$
   > (B^HB)^{-1} = \frac{1}{40-16}\begin{pmatrix} 20 & 4 \\ 4 & 2 \end{pmatrix} = \frac{1}{24}\begin{pmatrix} 20 & 4 \\ 4 & 2 \end{pmatrix}
   > $$
   >
   > $$
   > B^H = \begin{pmatrix} 1 & 0 & -1 \\ 0 & 2 & 4 \end{pmatrix},\quad C^H = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ -1 & 1 \\ 1 & 1 \end{pmatrix}
   > $$
   >
   > 先算
   > $$
   > (B^HB)^{-1}B^H = \frac{1}{24}\begin{pmatrix} 20 & 4 \\ 4 & 2 \end{pmatrix}\begin{pmatrix} 1 & 0 & -1 \\ 0 & 2 & 4 \end{pmatrix} = \frac{1}{24}\begin{pmatrix} 20 & 8 & -4 \\ 4 & 4 & 4 \end{pmatrix}
   > $$
   >
   > 再算
   > $$
   > C^H(CC^H)^{-1} = \frac{1}{3}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ -1 & 1 \\ 1 & 1 \end{pmatrix}
   > $$
   >
   > $$
   > A^+ = \frac{1}{72}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ -1 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} 20 & 8 & -4 \\ 4 & 4 & 4 \end{pmatrix} \\= \frac{1}{72}\begin{pmatrix} 20 & 8 & -4 \\ 4 & 4 & 4 \\ -16 & -4 & 8 \\ 24 & 12 & 0 \end{pmatrix} = \frac{1}{18}\begin{pmatrix} 5 & 2 & -1 \\ 1 & 1 & 1 \\ -4 & -1 & 2 \\ 6 & 3 & 0 \end{pmatrix}
   > $$
   >
   > **步骤三：判断相容性。**（用 $A^+$ 更方便，因为 $A^+ \in A\{1\}$）
   > $$
   > A^+b = \frac{1}{18}\begin{pmatrix} 5 & 2 & -1 \\ 1 & 1 & 1 \\ -4 & -1 & 2 \\ 6 & 3 & 0 \end{pmatrix} \begin{pmatrix} 4 \\ -2 \\ -2 \end{pmatrix} = \frac{1}{18}\begin{pmatrix} 20 - 4 + 2 \\ 4 - 2 - 2 \\ -16 + 2 - 4 \\ 24 - 6 + 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 0 \\ -1 \\ 1 \end{pmatrix}
   > $$
   >
   > $$
   > AA^+b = \begin{pmatrix} 1 & 0 & -1 & 1 \\ 0 & 2 & 2 & 2 \\ -1 & 4 & 5 & 3 \end{pmatrix} \begin{pmatrix} 1 \\ 0 \\ -1 \\ 1 \end{pmatrix} = \begin{pmatrix} 3 \\ 0 \\ -3 \end{pmatrix} \neq \begin{pmatrix} 4 \\ -2 \\ -2 \end{pmatrix} = b
   > $$
   >
   > $AA^+b \neq b$，故 **方程组不相容**。
   >
   > （也可用初等变换验证：增广矩阵秩为 $3 > \text{rank}(A) = 2$，不相容。）
   >
   > **步骤四：求极小最小二乘解。**
   >
   > 由定理1.6.5，对不相容方程组，$x = A^+b$ 是唯一的极小最小二乘解（唯一最佳逼近解）：
   > $$
   > x = A^+b = \begin{pmatrix} 1 \\ 0 \\ -1 \\ 1 \end{pmatrix}
   > $$
   >
   > 即是不相容方程组 $Ax = b$ 的极小最小二乘解。
   
9. 设 $A(t) = \begin{pmatrix} e^{2t} & t e^t \\ 0 & e^{-t} \end{pmatrix}$，
   （1）证明 $A(t)$ 可逆并求 $A^{-1}(t)$；

   （2）求 $\frac{\mathrm{d}}{\mathrm{d}t} A^{-1}(t)$。

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **命题2.1.1（求导法则）**（见《矩阵理论与应用II复习》 §1.2.1）：若 $A(t)$ 可逆，则 $\frac{d}{dt}A^{-1} = -A^{-1}\frac{dA}{dt}A^{-1}$。
   >
   > **（1）** $\det(A(t)) = e^{2t} \cdot e^{-t} - 0 \cdot (t e^t) = e^t \neq 0$（$\forall t$），故 $A(t)$ 可逆。
   >
   > 由二阶矩阵求逆公式 $\begin{pmatrix} a & b \\ c & d \end{pmatrix}^{-1} = \frac{1}{ad-bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$：
   > $$
   > A^{-1}(t) = \frac{1}{e^t}\begin{pmatrix} e^{-t} & -t e^t \\ 0 & e^{2t} \end{pmatrix} = \begin{pmatrix} e^{-2t} & -t \\ 0 & e^t \end{pmatrix}
   > $$
   >
   > **（2）** 利用命题2.1.1(5)：$\frac{d}{dt}A^{-1} = -A^{-1}\frac{dA}{dt}A^{-1}$。
   >
   > 先求 $\frac{dA}{dt} = \begin{pmatrix} 2e^{2t} & e^t + t e^t \\ 0 & -e^{-t} \end{pmatrix}$。
   >
   > $$
   > A^{-1}\frac{dA}{dt} = \begin{pmatrix} e^{-2t} & -t \\ 0 & e^t \end{pmatrix} \begin{pmatrix} 2e^{2t} & e^t + t e^t \\ 0 & -e^{-t} \end{pmatrix} \\
   > = \begin{pmatrix} 2 & e^{-t} + t e^{-t} + t e^{-t} \\ 0 & -1 \end{pmatrix} = \begin{pmatrix} 2 & e^{-t}(1+2t) \\ 0 & -1 \end{pmatrix}
   > $$
   >
   > $$
   > \frac{d}{dt}A^{-1} = -\begin{pmatrix} 2 & e^{-t}(1+2t) \\ 0 & -1 \end{pmatrix} \begin{pmatrix} e^{-2t} & -t \\ 0 & e^t \end{pmatrix} \\
   > = -\begin{pmatrix} 2e^{-2t} & -2t + e^{-t}(1+2t)e^t \\ 0 & -e^t \end{pmatrix}\\
   > = -\begin{pmatrix} 2e^{-2t} & -2t + 1 + 2t \\ 0 & -e^t \end{pmatrix} \\
   > = \begin{pmatrix} -2e^{-2t} & -1 \\ 0 & e^t \end{pmatrix}
   > $$
   >
   > 或直接对 $A^{-1}(t) = \begin{pmatrix} e^{-2t} & -t \\ 0 & e^t \end{pmatrix}$ 逐元素求导验证：
   >
   > $$
   > \frac{d}{dt}A^{-1} = \begin{pmatrix} -2e^{-2t} & -1 \\ 0 & e^t \end{pmatrix}
   > $$
   > 结果一致。

10. 设 $A(t) = \begin{pmatrix} t \cos t & e^t \sin t \\ t^2 + 1 & \ln(1+t) \end{pmatrix}$，计算：
    （1） $\frac{\mathrm{d}A(t)}{\mathrm{d}t}$；

    （2） $\int_0^t A(s) \mathrm{d}s$。

    > [!IMPORTANT]
    >
    > **【原理依据】** 使用 **定义2.1.1（函数矩阵的导数与积分）**（见《矩阵理论与应用II复习》 §1.2.1）：函数矩阵的导数和积分分别对每个元素逐项求导/积分即可。
    >
    > **（1）导数：**
    > $$
    > \frac{dA(t)}{dt} = \begin{pmatrix} \frac{d}{dt}(t \cos t) & \frac{d}{dt}(e^t \sin t) \\[4pt] \frac{d}{dt}(t^2+1) & \frac{d}{dt}(\ln(1+t)) \end{pmatrix}
    > $$
    >
    > - $\frac{d}{dt}(t \cos t) = \cos t - t \sin t$
    > - $\frac{d}{dt}(e^t \sin t) = e^t \sin t + e^t \cos t = e^t(\sin t + \cos t)$
    > - $\frac{d}{dt}(t^2+1) = 2t$
    > - $\frac{d}{dt}(\ln(1+t)) = \frac{1}{1+t}$
    >
    > $$
    > \frac{dA(t)}{dt} = \begin{pmatrix} \cos t - t \sin t & e^t(\sin t + \cos t) \\ 2t & \frac{1}{1+t} \end{pmatrix}
    > $$
    >
    > **（2）积分：** 
    > $$
    > \int_0^t A(s)ds = \begin{pmatrix} \int_0^t s \cos s \, ds & \int_0^t e^s \sin s \, ds \\[4pt] \int_0^t (s^2+1)ds & \int_0^t \ln(1+s)ds \end{pmatrix}
    > $$
    >
    > - $\int_0^t s \cos s \, ds = [s \sin s]_0^t - \int_0^t \sin s \, ds = t \sin t + \cos t - 1$
    > - $\int_0^t e^s \sin s \, ds = \frac{e^s}{2}(\sin s - \cos s)\big|_0^t = \frac{e^t}{2}(\sin t - \cos t) + \frac{1}{2}$
    > - $\int_0^t (s^2+1)ds = \frac{t^3}{3} + t$
    > - $\int_0^t \ln(1+s)ds = [(1+s)\ln(1+s) - s]_0^t = (1+t)\ln(1+t) - t$
    >
    > $$
    > \int_0^t A(s)ds = \begin{pmatrix} t\sin t + \cos t - 1 & \frac{e^t}{2}(\sin t - \cos t) + \frac{1}{2} \\[4pt] \frac{t^3}{3} + t & (1+t)\ln(1+t) - t \end{pmatrix}
    > $$
    
11. 设 $A = \begin{pmatrix} -1 & -2 & 6 \\ -1 & 0 & 3 \\ -1 & -1 & 4 \end{pmatrix}, \mathbf{b}(t) = \begin{pmatrix} e^t \\ 2e^t \\ e^t \end{pmatrix}$，求 $e^{At}$，并求微分方程 $\frac{\mathrm{d}\mathbf{x}(t)}{\mathrm{d}t} = A\mathbf{x}(t) + \mathbf{b}(t)$ 满足初始条件 $\mathbf{x}(0) = (-1, 0, 1)^{\mathrm{T}}$ 的解。

    > [!IMPORTANT]
    >
    > **【原理依据】** 使用 **谱上一致性方法**（见《矩阵论复习》 第35条 & 第39条）求 $e^{At}$，再使用一阶线性常系数微分方程组的通解公式，即常数变异法（见《矩阵理论与应用II复习》 §1.2.3 应用二）：
    > $$
    > x(t) = e^{At}x(0) + \int_0^t e^{A(t-s)}b(s)ds
    > $$
    >
    > ---
    >
    > **第一部分：求 $e^{At}$（谱上一致性方法）**
    >
    > **步骤一：** 求 $A$ 的特征多项式。
    > $$
    > |\lambda I - A| = \begin{vmatrix} \lambda+1 & 2 & -6 \\ 1 & \lambda & -3 \\ 1 & 1 & \lambda-4 \end{vmatrix}
    > $$
    > 按第一行展开：
    > $$
    > = (\lambda+1)\begin{vmatrix} \lambda & -3 \\ 1 & \lambda-4 \end{vmatrix} - 2\begin{vmatrix} 1 & -3 \\ 1 & \lambda-4 \end{vmatrix} + (-6)\begin{vmatrix} 1 & \lambda \\ 1 & 1 \end{vmatrix}\\
    > =(\lambda+1)[\lambda(\lambda-4) + 3] - 2[(\lambda-4) + 3] - 6[1 - \lambda]\\
    > = (\lambda+1)(\lambda^2 - 4\lambda + 3) - 2(\lambda - 1) - 6(1 - \lambda)\\
    > = (\lambda+1)(\lambda-1)(\lambda-3) - 2(\lambda-1) + 6(\lambda-1)\\
    > = (\lambda-1)[(\lambda+1)(\lambda-3) + 4] = (\lambda-1)[\lambda^2 - 2\lambda + 1] = (\lambda-1)^3=0
    > $$
    > 故 $A$ 的特征值为 $\lambda = 1$（三重根）。
    >
    > **步骤二：** 求最小多项式。计算 $(A - I)^2$：
    >
    > $$
    > A - I = \begin{pmatrix} -2 & -2 & 6 \\ -1 & -1 & 3 \\ -1 & -1 & 3 \end{pmatrix}
    > $$
    >
    > $$
    > (A - I)^2 = \begin{pmatrix} -2 & -2 & 6 \\ -1 & -1 & 3 \\ -1 & -1 & 3 \end{pmatrix}^2 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix} = O
    > $$
    >
    > 故最小多项式 $m_A(\lambda) = (\lambda - 1)^2$（二次，因为 $A - I \neq O$ 但 $(A-I)^2 = O$）。
    >
    > **步骤三：** 设余式多项式。因 $m_A(\lambda)$ 为二次，设 $r(\lambda) = \alpha_0(t) + \alpha_1(t)\lambda$。
    >
    > 令 $f(\lambda) = e^{\lambda t} = m_A(\lambda)g(\lambda) + r(\lambda)$。在 $\lambda = 1$（二重根）处，有：
    > $$
    > \begin{cases} r(1) = f(t) \Rightarrow \alpha_0 + \alpha_1 = e^t \\ r'(1) = f'(t) \Rightarrow \alpha_1 = t e^t \end{cases}
    > $$
    >
    > 解得 $\alpha_1 = t e^t$，$\alpha_0 = e^t - t e^t = (1-t)e^t$。
    >
    > **步骤四：**
    > $$
    > e^{At} = r(A) = \alpha_0 I + \alpha_1 A = (1-t)e^t I + t e^t A = e^t[(1-t)I + tA]
    > $$
    >
    > $$
    > = e^t \left[ \begin{pmatrix} 1-t & 0 & 0 \\ 0 & 1-t & 0 \\ 0 & 0 & 1-t \end{pmatrix} + \begin{pmatrix} -t & -2t & 6t \\ -t & 0 & 3t \\ -t & -t & 4t \end{pmatrix} \right]
    > $$
    >
    > $$
    > e^{At} = e^t \begin{pmatrix} 1-2t & -2t & 6t \\ -t & 1-t & 3t \\ -t & -t & 1+3t \end{pmatrix}
    > $$
    >
    > ---
    >
    > **第二部分：求微分方程的解**
    >
    > 由通解公式 $x(t) = e^{At}x(0) + \int_0^t e^{A(t-s)}b(s)ds$。
    >
    > **齐次部分：**
    >
    > $$
    > e^{At}x(0) = e^t \begin{pmatrix} 1-2t & -2t & 6t \\ -t & 1-t & 3t \\ -t & -t & 1+3t \end{pmatrix} \begin{pmatrix} -1 \\ 0 \\ 1 \end{pmatrix} \\= e^t \begin{pmatrix} -(1-2t) + 6t \\ t + 3t \\ t + (1+3t) \end{pmatrix} \\= e^t \begin{pmatrix} -1 + 8t \\ 4t \\ 1 + 4t \end{pmatrix}
    > $$
    >
    > **非齐次部分（特解积分）：** 先计算 $e^{A(t-s)}b(s)$。
    >
    > $$
    > e^{A(t-s)} = e^{t-s} \begin{pmatrix} 1-2(t-s) & -2(t-s) & 6(t-s) \\ -(t-s) & 1-(t-s) & 3(t-s) \\ -(t-s) & -(t-s) & 1+3(t-s) \end{pmatrix}
    > $$
    >
    > $$
    > b(s) = \begin{pmatrix} e^s \\ 2e^s \\ e^s \end{pmatrix} = e^s \begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}
    > $$
    >
    > 注意到 $e^{A(t-s)}b(s) = e^t \cdot M(t-s) \begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}$，其中 $M(\tau) = \begin{pmatrix} 1-2\tau & -2\tau & 6\tau \\ -\tau & 1-\tau & 3\tau \\ -\tau & -\tau & 1+3\tau \end{pmatrix}$。
    >
    > 计算 $M(t-s)\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}$：
    > $$
    > \begin{pmatrix} (1-2\tau) \cdot 1 + (-2\tau) \cdot 2 + 6\tau \cdot 1 \\ (-\tau) \cdot 1 + (1-\tau) \cdot 2 + 3\tau \cdot 1 \\ (-\tau) \cdot 1 + (-\tau) \cdot 2 + (1+3\tau) \cdot 1 \end{pmatrix} = \begin{pmatrix} 1 - 2\tau - 4\tau + 6\tau \\ -\tau + 2 - 2\tau + 3\tau \\ -\tau - 2\tau + 1 + 3\tau \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}
    > $$
    >
    > 惊喜！$M(\tau)\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}$（$\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}$ 恰是 $A$ 的特征向量，对应特征值 $1$）。故：
    >
    > $$
    > e^{A(t-s)}b(s) = e^t \begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} = b(t)
    > $$
    >
    > $$
    > \int_0^t e^{A(t-s)}b(s)ds = \int_0^t b(t) ds = t \cdot b(t) = t e^t \begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix}
    > $$
    >
    > **最终解：**
    > $$
    > x(t) = e^{At}x(0) + \int_0^t e^{A(t-s)}b(s)ds \\= e^t\begin{pmatrix} -1 + 8t \\ 4t \\ 1 + 4t \end{pmatrix} + t e^t\begin{pmatrix} 1 \\ 2 \\ 1 \end{pmatrix} \\= e^t\begin{pmatrix} -1 + 9t \\ 6t \\ 1 + 5t \end{pmatrix}
    > $$
    >
    > 验证初始条件：$x(0) = e^0\begin{pmatrix} -1 \\ 0 \\ 1 \end{pmatrix} = \begin{pmatrix} -1 \\ 0 \\ 1 \end{pmatrix}$。

12. 证明：若 $A\mathbf{X} = \mathbf{b}$ 不相容，则 $A^H A \mathbf{X} = A^H \mathbf{b}$ 一定有解，且它的解也是 $A\mathbf{X} = \mathbf{b}$ 的最小二乘解。

    > [!IMPORTANT]
    >
    > **【原理依据】** 使用 **推论1.5.3** 与 **注2**（见《矩阵理论与应用II复习》 §1.1.5）：方程组 $A^HAx = A^Hb$ 与 $Ax = AA_l^-b$ 是相容同解方程组，且其解即为不相容方程组 $Ax = b$ 的最小二乘解。
    >
    > **【证明】**
    >
    > **（1）$A^HAx = A^Hb$ 一定有解（相容性）：**
    >
    > $A^Hb \in R(A^H)$，而 $R(A^HA) = R(A^H)$（由引理1.1.1，$\text{rank}(A^HA) = \text{rank}(A) = \text{rank}(A^H)$，且 $R(A^HA) \subseteq R(A^H)$，维数相同故相等）。因此 $A^Hb \in R(A^HA)$，即存在 $x$ 使得 $A^HAx = A^Hb$，方程组相容。
    >
    > **（2）$A^HAx = A^Hb$ 的解是 $Ax = b$ 的最小二乘解：**
    >
    > 最小二乘问题 $\min_x\|Ax - b\|_2$ 的极值条件（见《矩阵理论与应用II复习》 §1.2.3 应用一）为：
    > $$
    > \frac{d}{dx}\|Ax - b\|_2^2 = 2A^HAx - 2A^Hb = 0 \Rightarrow A^HAx = A^Hb
    > $$
    >
    > 此即正规方程。由正交投影的角度（见《矩阵理论与应用II复习》 §1.1.5 推论1.5.3 及注2），$Ax = b$ 的最小二乘解满足 $Ax = AA_l^-b$（$b$ 在 $R(A)$ 上的正交投影），而 $A^HAx = A^Hb$ 等价于 $A^H(Ax - b) = 0$，即残量 $Ax - b \perp R(A)$，这正是最佳逼近的充要条件（**定理3.1.4**，见《矩阵理论与应用II复习》 §1.3.2）。
    >
    > 故 $A^HAx = A^Hb$ 的解确为 $Ax = b$ 的最小二乘解。
    >
    > **注：** 由《矩阵理论与应用II复习》 §1.1.5 注2，方程组 $A^HAx = A^Hb$、$Ax = AA_l^-b$ 和 $Ax = A(A^HA)^-A^Hb$ 是相容同解方程组，它们都等价于最小二乘问题。

# 2.往年真题

这一部分只整理一些不太常规的出题方式、易错点以及一些比较细节的结论性知识点。                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

1. 设$A$是$n$阶可逆矩阵，则$\begin{pmatrix} A & A \\ A & A \end{pmatrix}$的伪逆是

   > [!IMPORTANT]
   >
   > **【原理依据】** 使用 **定理1.6.4（满秩分解求加号逆）**（见《矩阵理论与应用II复习》 §1.1.6）：若 $B = CD$ 为满秩分解，则 $B^+ = D^H(DD^H)^{-1}(C^HC)^{-1}C^H$。再利用加号逆唯一性（**定理1.6.2**），通过验证四条Penrose方程确认结果。
   >
   > ---
   >
   > **【解法一：满秩分解法（推荐）】**
   >
   > 记 $B = \begin{pmatrix} A & A \\ A & A \end{pmatrix} \in \mathbb{C}^{2n \times 2n}$。由于 $A$ 可逆，易见：
   >
   > $$
   > B = \begin{pmatrix} A \\ A \end{pmatrix} \begin{pmatrix} I_n & I_n \end{pmatrix}
   > $$
   >
   > 令 $C = \begin{pmatrix} A \\ A \end{pmatrix} \in \mathbb{C}^{2n \times n}$（列满秩，因 $A$ 可逆），$D = \begin{pmatrix} I_n & I_n \end{pmatrix} \in \mathbb{C}^{n \times 2n}$（行满秩），则 $B = CD$ 是 $B$ 的满秩分解，$\text{rank}(B) = n$。
   >
   > **步骤一：计算 $D^H(DD^H)^{-1}$**
   >
   > $$
   > D^H = \begin{pmatrix} I_n \\ I_n \end{pmatrix}, DD^H = \begin{pmatrix} I_n & I_n \end{pmatrix}\begin{pmatrix} I_n \\ I_n \end{pmatrix} = 2I_n,  (DD^H)^{-1} = \frac{1}{2}I_n
   > $$
   >
   > 
   > $$
   > D^H(DD^H)^{-1} = \frac{1}{2}\begin{pmatrix} I_n \\ I_n \end{pmatrix}
   > $$
   >
   > **步骤二：计算 $(C^HC)^{-1}C^H$**
   >$$
   > C^H = \begin{pmatrix} A^H & A^H \end{pmatrix}, C^HC = \begin{pmatrix} A^H & A^H \end{pmatrix}\begin{pmatrix} A \\ A \end{pmatrix} = A^HA + A^HA = 2A^HA
   > $$
   > 
   >$$
   > (C^HC)^{-1} = \frac{1}{2}(A^HA)^{-1} = \frac{1}{2}A^{-1}(A^H)^{-1}
   > $$
   > 
   >$$
   > (C^HC)^{-1}C^H = \frac{1}{2}A^{-1}(A^H)^{-1}\begin{pmatrix} A^H & A^H \end{pmatrix} = \frac{1}{2}A^{-1}\begin{pmatrix} I_n & I_n \end{pmatrix}
   > $$
   > 
   >**步骤三：合并求 $B^+$**
   > 
   >$$
   > B^+ = D^H(DD^H)^{-1}(C^HC)^{-1}C^H = \frac{1}{2}\begin{pmatrix} I_n \\ I_n \end{pmatrix} \cdot \frac{1}{2}A^{-1}\begin{pmatrix} I_n & I_n \end{pmatrix}
   > $$
   > 
   >$$
   > B^+ = \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix}
   > $$
   > 
   >---
   > 
   >**【解法二：直接验证四条Penrose方程】**
   > 
   >猜想 $B^+ = \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix}$，逐条验证：
   > 
   >**（1）$BB^+B = B$：**
   > 
   >$$
   > BB^+ = \begin{pmatrix} A & A \\ A & A \end{pmatrix} \cdot \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix} = \frac{1}{4}\begin{pmatrix} 2I_n & 2I_n \\ 2I_n & 2I_n \end{pmatrix} = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix}
   > $$
   > 
   >$$
   > BB^+B = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix}\begin{pmatrix} A & A \\ A & A \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 2A & 2A \\ 2A & 2A \end{pmatrix} = \begin{pmatrix} A & A \\ A & A \end{pmatrix} = B
   > $$
   > 
   >满足第一条 Penrose 方程。
   > 
   >**（2）$B^+BB^+ = B^+$：**
   > 
   >$$
   > B^+B = \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix}\begin{pmatrix} A & A \\ A & A \end{pmatrix} = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix}
   > $$
   > 
   >$$
   > B^+BB^+ = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix} \cdot \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix} = \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix} = B^+
   > $$
   > 
   >满足第二条 Penrose 方程。
   > 
   >**（3）$(BB^+)^H = BB^+$：**
   > 
   >$$
   > BB^+ = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix}, \quad (BB^+)^H = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix} = BB^+
   > $$
   > 
   >满足第三条 Penrose 方程。
   > 
   >**（4）$(B^+B)^H = B^+B$：**
   > 
   >$B^+B = \frac{1}{2}\begin{pmatrix} I_n & I_n \\ I_n & I_n \end{pmatrix}$ 显然为 Hermite 矩阵，$(B^+B)^H = B^+B$。满足第四条 Penrose 方程。
   > 
   >由加号逆唯一性（**定理1.6.2**），
   > 
   >$$
   > B^+ = \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix}
   > $$
   > 
   >---
   > 
   >**【注】** 此题也可利用推论1.6.2：$C$ 列满秩时 $C^+ = (C^HC)^{-1}C^H$，$D$ 行满秩时 $D^+ = D^H(DD^H)^{-1}$，结合满秩分解 $B=CD$ 直接求得 $B^+ = D^+C^+$（注意：一般情况下 $(CD)^+ \neq D^+C^+$，但此处因 $C$ 列满秩、$D$ 行满秩且乘积结构与秩匹配，可验证等号成立）。实际上此题中 $D^+C^+ = \frac{1}{2}\begin{pmatrix} I_n \\ I_n \end{pmatrix} \cdot \frac{1}{2}A^{-1}\begin{pmatrix} I_n & I_n \end{pmatrix} = \frac{1}{4}\begin{pmatrix} A^{-1} & A^{-1} \\ A^{-1} & A^{-1} \end{pmatrix}$，结果一致。

2. 设$A=\begin{pmatrix} 2 & i \\ 0 & 0  \\ 0 & 0 \end{pmatrix}$，求$A$的奇异值分解。

   > [!WARNING]
   >
   > 注意，对于带有复数的矩阵，转置时需要将每一项中的虚部取相反数。
   >
   > 本题中，$A^H=\begin{pmatrix} 2 & 0 & 0 \\ -i & 0 & 0 \end{pmatrix}$。一定注意，不然后续会计算错误。

3. $A=\begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$，$B=\begin{pmatrix} a & 1 & 2 \\ 0 & b & 1 \\ 0&0&c\end{pmatrix}$，则$A\otimes B=$__________________，$A\otimes B$的特征值为__________________。

   > [!IMPORTANT]
   >
   > 第一空使用矩阵直积的定义，第二空，矩阵直积的特征值为

4. 已知$A$的广义逆$A^+=\begin{pmatrix} 1 & 3 & -2 \\ -\frac{3}{2} & -3 & \frac{5}{2} \\ 1&1&-1\end{pmatrix}$，求矩阵$A$。

   > [!IMPORTANT]
   >
   > **【原理依据】**
   > - **命题1.6.1(1)**（见《矩阵理论与应用II复习》 §1.1.6）：$(A^+)^+ = A$；
   > - **引理**：若 $n$ 阶方阵 $G$ 可逆，则 $G^+ = G^{-1}$（证明见下文步骤二）；
   > - **定理1.6.2（加号逆唯一性）**：满足四条 Penrose 方程的矩阵唯一。
   >
   > **【逻辑链条】**
   > $$A^+ \text{ 可逆} \;\xrightarrow{\text{引理}}\; (A^+)^+ = (A^+)^{-1} \;\xrightarrow{\text{命题1.6.1(1)}}\; A = (A^+)^+ = (A^+)^{-1}$$
   >
   > 故只需验证 $A^+$ 可逆，再求其普通逆即可。
   >
   > ---
   >
   > **步骤一：验证 $A^+$ 可逆。**
   >
   > 计算行列式：
   >
   > $$
   > \begin{aligned}
   > \det(A^+) &= \begin{vmatrix} 1 & 3 & -2 \\[4pt] -\frac{3}{2} & -3 & \frac{5}{2} \\[4pt] 1 & 1 & -1 \end{vmatrix} \\[8pt]
   > &= 1 \cdot \begin{vmatrix} -3 & \frac{5}{2} \\[4pt] 1 & -1 \end{vmatrix} - 3 \cdot \begin{vmatrix} -\frac{3}{2} & \frac{5}{2} \\[4pt] 1 & -1 \end{vmatrix} + (-2) \cdot \begin{vmatrix} -\frac{3}{2} & -3 \\[4pt] 1 & 1 \end{vmatrix} \\[8pt]
   > &= 1 \cdot \left(3 - \frac{5}{2}\right) - 3 \cdot \left(\frac{3}{2} - \frac{5}{2}\right) - 2 \cdot \left(-\frac{3}{2} + 3\right) \\[8pt]
   > &= \frac{1}{2} - 3 \cdot (-1) - 2 \cdot \frac{3}{2} \\[8pt]
   > &= \frac{1}{2} + 3 - 3 = \frac{1}{2} \neq 0
   > \end{aligned}
   > $$
   >
   > 故 $A^+$ 可逆，$\text{rank}(A^+) = 3$。
   >
   > ---
   >
   > **步骤二：证明引理——若 $n$ 阶方阵 $G$ 可逆，则 $G^+ = G^{-1}$。**
   >
   > 由 **定理1.6.2（加号逆唯一性定理）**，加号逆是唯一的。因此只需验证 $G^{-1}$ 满足四条 Penrose 方程即可：
   >
   > | Penrose 方程 | 验证 |
   > |:---:|:---|
   > | (1) $G X G = G$ | $G \cdot G^{-1} \cdot G = I \cdot G = G$ ✓ |
   > | (2) $X G X = X$ | $G^{-1} \cdot G \cdot G^{-1} = I \cdot G^{-1} = G^{-1}$ ✓ |
   > | (3) $(GX)^H = GX$ | $(G G^{-1})^H = I^H = I = G G^{-1}$ ✓ |
   > | (4) $(XG)^H = XG$ | $(G^{-1} G)^H = I^H = I = G^{-1} G$ ✓ |
   >
   > $G^{-1}$ 全部满足，由唯一性即得 $G^+ = G^{-1}$。
   >
   > **应用于本题：** 令 $G = A^+$，则 $(A^+)^+ = (A^+)^{-1}$。再由命题1.6.1(1)：$(A^+)^+ = A$，故 $A = (A^+)^{-1}$。
   >
   > ---
   >
   > **步骤三：求 $(A^+)^{-1}$。**
   >
   > 使用伴随矩阵法。先计算各元素的代数余子式：
   >
   > $$
   > \begin{aligned}
   > C_{11} &= +\begin{vmatrix} -3 & \frac{5}{2} \\[4pt] 1 & -1 \end{vmatrix} = 3 - \frac{5}{2} = \frac{1}{2}, &
   > C_{12} &= -\begin{vmatrix} -\frac{3}{2} & \frac{5}{2} \\[4pt] 1 & -1 \end{vmatrix} = -\left( \frac{3}{2} - \frac{5}{2} \right) = 1, \\[8pt]
   > C_{13} &= +\begin{vmatrix} -\frac{3}{2} & -3 \\[4pt] 1 & 1 \end{vmatrix} = -\frac{3}{2} + 3 = \frac{3}{2}, &
   > C_{21} &= -\begin{vmatrix} 3 & -2 \\[4pt] 1 & -1 \end{vmatrix} = -(-3 + 2) = 1, \\[8pt]
   > C_{22} &= +\begin{vmatrix} 1 & -2 \\[4pt] 1 & -1 \end{vmatrix} = -1 + 2 = 1, &
   > C_{23} &= -\begin{vmatrix} 1 & 3 \\[4pt] 1 & 1 \end{vmatrix} = -(1 - 3) = 2, \\[8pt]
   > C_{31} &= +\begin{vmatrix} 3 & -2 \\[4pt] -3 & \frac{5}{2} \end{vmatrix} = \frac{15}{2} - 6 = \frac{3}{2}, &
   > C_{32} &= -\begin{vmatrix} 1 & -2 \\[4pt] -\frac{3}{2} & \frac{5}{2} \end{vmatrix} = -\left( \frac{5}{2} - 3 \right) = \frac{1}{2}, \\[8pt]
   > C_{33} &= +\begin{vmatrix} 1 & 3 \\[4pt] -\frac{3}{2} & -3 \end{vmatrix} = -3 + \frac{9}{2} = \frac{3}{2}.
   > \end{aligned}
   > $$
   >
   > 余子式矩阵 $\begin{pmatrix} \frac{1}{2} & 1 & \frac{3}{2} \\[4pt] 1 & 1 & 2 \\[4pt] \frac{3}{2} & \frac{1}{2} & \frac{3}{2} \end{pmatrix}$，转置得伴随矩阵：
   >
   > $$
   > \text{adj}(A^+) = \begin{pmatrix} \frac{1}{2} & 1 & \frac{3}{2} \\[4pt] 1 & 1 & \frac{1}{2} \\[4pt] \frac{3}{2} & 2 & \frac{3}{2} \end{pmatrix}
   > $$
   >
   > 故：
   >
   > $$
   > A = (A^+)^{-1} = \frac{1}{\det(A^+)} \cdot \text{adj}(A^+) = 2 \cdot \begin{pmatrix} \frac{1}{2} & 1 & \frac{3}{2} \\[4pt] 1 & 1 & \frac{1}{2} \\[4pt] \frac{3}{2} & 2 & \frac{3}{2} \end{pmatrix}
   > $$
   >
   > $$
   > A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 2 & 1 \\ 3 & 4 & 3 \end{pmatrix}
   > $$
   >
   > ---
   >
   > **步骤四：验证 $A^+ A = I_3$。**
   >$$
   > \begin{aligned}
   > A^+A &= \begin{pmatrix} 1 & 3 & -2 \\[4pt] -\frac{3}{2} & -3 & \frac{5}{2} \\[4pt] 1 & 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 2 & 3 \\ 2 & 2 & 1 \\ 3 & 4 & 3 \end{pmatrix} \\[8pt]
   > &= \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} = I_3
   > \end{aligned}
   > $$
   > 
   > 验证通过。由 $A^+A = I_3$ 及 $A^+$ 可逆亦得 $AA^+ = I_3$，故求得的 $A$ 确实以给定矩阵为其加号逆。
   > 
   > $$
   >A = \begin{pmatrix} 1 & 2 & 3 \\ 2 & 2 & 1 \\ 3 & 4 & 3 \end{pmatrix}
   > $$

5. 













