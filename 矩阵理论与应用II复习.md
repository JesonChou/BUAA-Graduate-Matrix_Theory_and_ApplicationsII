<h1 align = "center">矩阵理论与应用II复习</h1>
<p align="right">SY2502110 朱宛瑜</p>

# 1.基本概念

## 1.1 广义逆矩阵及其应用

### 1.1.1 基本概念及奇异值分解

1. 定义1.1.1（广义逆矩阵）：设$A\in \mathbb{C}^{m×n}$，若$X\in \mathbb{C}^{n×m}$满足以下四个Penrose方程的全部或一部分：

   > [!IMPORTANT]
   >
   > 1、$AXA=A$
   >
   > 2、$XAX=X$
   >
   > 3、$(AX)^H=AX$
   >
   > 4、$(XA)^H=XA$

   则称矩阵$X$是矩阵$A$的==广义逆矩阵==。

   符号表示为$G=A^{(i,j,k)}$，其中$i,j,k=1,2,3,4$且$i,j,k$互不相等。

   满足全部四条的记为$G=A^{(1,2,3,4)}$或者$A^+$，称为加号逆、伪逆或者$Moore-Penrose$广义逆。

   $A$为奇异矩阵或者非方阵时，$A$的减号逆矩阵不一定唯一，记作集合$A\{1\}$。其他逆同样。

2. 常见广义逆为$A^{(1)}$，$A^{(1,3)}$，$A^{(1,4)}$和$A^+$。其中$A^{(1)}$称为减号逆，记作$A^-$；$A^{(1,3)}$称为最小二乘广义逆，记作$A^-_l$；$A^{(1,4)}$称为极小范数广义逆，记为$A_m^-$。

3. $Jordan$标准形分解：$J=P^-AP=diag(J_1,J_2,...,J_s)$。

4. $Jordan$标准形分解的局限性：

   > [!IMPORTANT]
   >
   > 1、矩阵$A$必须是方阵；
   >
   > 2、$Jordan$阵不一定是对角阵；
   >
   > 3、相似变换矩阵不是酉矩阵。
   >
   > 注：酉矩阵的定义：$A^HA=AA^H=I$。即复空间上的正交矩阵（$A^TA=AA^T=I$）。

   但是对于$A\in\mathbb{C}^{m×n}$，$m$不一定等于$n$，即$A$是==非方阵==。能否找到一种可以克服$Jordan$分解的局限性的分解方法？

5. 基于上述思路，我们选择Hermite矩阵（即$A^H=A$）$A^HA_{n×n}$与$AA^H_{m×m}$，则分别存在酉矩阵$U\in\mathbb{C}^{m×m}$与$V\in\mathbb{C}^{n×n}$，使得
   $$
   U^HAA^HU=\Sigma_1,V^HA^HAV=\Sigma_2
   $$
   式中$\Sigma_1$与$\Sigma_2$为对角矩阵，对角线元素为矩阵$AA^H$和$A^HA$的==特征值==。<font color='red'>（这里涉及到Hermite矩阵的性质，Hermite矩阵是正规矩阵，其一定可以被酉相似对角化）</font>

   则矩阵$A$可以分解为
   $$
   A=UDV^H
   $$
   其中$DD^H=\Sigma_1$，$D^HD=\Sigma_2$。

   > [!IMPORTANT]
   >
   > 具体推导步骤为：
   > $$
   > AA^H=UDV^H·(UDV^H)^H=UDV^HVD^HU^H
   > $$
   > 由于$U$和$V$均为酉矩阵，则根据酉矩阵的性质（$A^HA=AA^H=I$）可得
   > $$
   > AA^H=UDD^HU^H
   > $$
   > 将式$(4)$中的矩阵左乘$U^H$并右乘$U$（因为$U$是酉矩阵，所以$U^HU=UU^H=I$）：
   > $$
   > U^HAA^HU=DD^H=\Sigma_1
   > $$
   > 对于$V$和$A^HA$同理。

6. ==引理1.1.1== 设$A\in\mathbb{C}^{m×n}$，则$rank(A^HA)=rank(AA^H)=rank(A)$。

   > [!NOTE]
   >
   > （以下推导过程需要了解掌握，可能出证明题）
   >
   > 首先需要知道矩阵秩的性质，即$rank(A)=rank(A^H)$。
   >
   > 设$\vec x\in N(A^HA)$，则有$A^HA\vec{x}=0$。

7. ==引理1.1.2== 设$A\in\mathbb{C}^{m×n}$，则$A^HA$与$AA^H$都是半正定Hermite矩阵。

   > [!NOTE]
   >
   > 矩阵$A$是半正定Hermite矩阵还有以下命题等价（PPT P13）
   >
   
8. ==引理1.1.3== 设$A\in\mathbb{C}^{m×n}_r$，则$A^HA$与$AA^H$的所有非零特征值完全相同且**非零特征值的个数**均为$r$。<font color='red'>（这里与奇异值的定义相关，$A^HA$与$AA^H$的非零特征值个数为$r$，则矩阵$A$的正奇异值的个数为$r$）</font>

9. 计算$A^HA$和$AA^H$的特征值时，可优先考虑计算阶数较小的矩阵。

10. ==定义1.1.2（奇异值）==设$A\in\mathbb{C}_r^{m×n}$，$A^HA$的特征值满足$\lambda_1\geq\lambda_2\geq···\geq\lambda_r>0$，$\lambda_{r+1}=\lambda_{r+2}=···=\lambda_n=0$，称$\sigma_i=\sqrt{\lambda_i}$为$A$的奇异值。称$\sigma_i,i=1,···,r$为$A$的正奇异值。

11. 由引理1.1.3可知，要计算一个矩阵$A$的正奇异值，只需要计算$A^HA$或者$AA^H$即可。

12. ==思考==：复方阵的特征值和奇异值有何异同？

    同一矩阵的特征值可能与其奇异值完全相同，可能部分相同，也可能完全不同。

13. ==命题1.1.1== 设$\sigma_1\geq\sigma_2\geq···\geq \sigma_n$和$|\lambda_1|\geq|\lambda_2|\geq···\geq|\lambda_n|$分别是$n$阶正规矩阵$A$（$A^HA=AA^H$）的奇异值和特征值，则$\sigma_i=|\lambda_i|,i=1,···,n$。也就是说正规矩阵的奇异值依顺序等于其特征值的绝对值。

14. ==命题1.1.2== 设$\sigma_i$和$\lambda_i(i=1,···,n)$分别是$n$阶复方阵$A$的奇异值和特征值，则
    $$
    \sum_{i=1}^n|\lambda_i|^2\leq\sum_{i=1}^n|\sigma_i|^2
    $$
    该命题的证明借助了Schur引理。

15. ==定理1.1.1（奇异值分解）==设$A\in\mathbb{C}_r^{m×n}$，则存在酉矩阵$U\in\mathbb{C}^{m×m}$与酉矩阵$V\in\mathbb{C}^{n×n}$使得
    $$
    A=U\begin{bmatrix}
    \Sigma_r & 0 \\
    0 & 0
    \end{bmatrix}V^H=UDV^H
    $$
    其中$\Sigma_r=diag(\sigma_1,···,\sigma_r)$，$\sigma_i(i=1,···,r)$是矩阵$A$的正奇异值。

    > [!CAUTION]
    >
    > 若$A=UDV^H$是奇异值分解，则酉矩阵$U$的列向量为$AA^H$的特征向量，酉矩阵$V$的列向量为$A^HA$的特征向量。

16. 奇异值分解步骤1：

    > [!IMPORTANT]
    >
    > 1、计算$AA^H$的$m$个奇异值，确定对角阵$\Sigma_r$；
    >
    > 2、计算$AA^H$的$m$个标准正交特征向量$\alpha_1,···,\alpha_m$，将其标准化构成酉矩阵$U$；（这一步需要采用诸如施密特正交化的方法进行）
    >
    > 3、计算$A^H\alpha_1,···,A^H\alpha_r$，将其扩为$\mathbb{C}^n$的一组标准正交基，构成矩阵$V$。
    >
    > $A^H\alpha_1,···,A^H\alpha_r$实际上是$A^HA$对应于$\lambda_1,···,\lambda_r$的正交特征向量，其推导步骤如下：
    >
    > $A^HA(A^H\alpha_1)=A^H(AA^H\alpha_1)=A^H(\lambda_1\alpha_1)=\lambda_1(A^H\alpha_1)$。满足特征向量的定义。

17. 奇异值分解步骤2：

    > [!IMPORTANT]
    >
    > 1、计算$A^HA$的$n$个奇异值，确定对角阵$\Sigma_r$；
    >
    > 2、计算$A^HA$的$n$个标准正交特征向量$\beta_1,···,\beta_n$，将其标准话构成酉矩阵$V$；
    >
    > 3、计算$A\beta_1/\sigma_1,···,A\beta_r/\sigma_r$，将其扩展为$\mathbb{C}^m$的一组标准正交基，构成矩阵$U$。
    >
    > $A\beta_1/\sigma_1,···,A\beta_r/\sigma_r$实际上是$AA^H$对应于$\lambda_1,···,\lambda_r$的正交特征向量。

18. ==扩充标准正交基==：在求解奇异值分解时，对于酉矩阵$U$或者$V$，有时会出现三维空间只有两个正交基向量。此时需要将其扩充为$\mathbb{R}^3$中的标准正交基，用以满足酉矩阵的定义。

19. 尽管矩阵$U$的列向量为矩阵$AA^H$的特征向量，矩阵$V$的列向量为$A^HA$的特征向量，但是酉矩阵$U$和$V$的选择具有一定相关性。可以先确定酉矩阵$U$，再根据$V_1=A^HU_1\Sigma_r^{-1}$确定$V$的前$r$列；或者先确定酉矩阵$V$，根据$U_1=AV_1\Sigma_r^{-1}$确定$U$的前$r$列。

20. ==定理1.1.2（极分解）==任意复方阵$A$必有如下分解：
    $$
    A=GW
    $$
    式中，$G$为半正定Hermite矩阵，$W$是酉矩阵。当矩阵$A$可逆时，$G$是正定Hermite矩阵，此时该极分解式唯一。

    上述分解式称为左极分解。若$A=VG$，其中$V$是酉矩阵，$G$是半正定Hermite矩阵，则称为右极分解。

21. ==应用（最小二乘问题）==：奇异值分解在矩阵特征值、广义逆矩阵等矩阵分析和计算方面有着重要应用，而且在图像处理、机器学习等领域有着广泛应用。

### 1.1.2 矩阵方程$AXB=D$

1. ==定义1.2.1（相容矩阵）==已知矩阵$A\in\mathbb{C}_{r_A}^{m×n}$，$B\in\mathbb{C}_{r_B}^{p×q}$，$D\in\mathbb{C}_{r_D}^{m×q}$，待定矩阵$X\in\mathbb{C}_{r_X}^{n×p}$，定义矩阵方程
   $$
   AXB=D
   $$
   若存在矩阵$X\in\mathbb{C}_{r_X}^{n×p}$使矩阵方程$AXB=D$成立，则称方程$AXB=D$**相容**，否则称为**不相容方程**或者**矛盾方程**。

2. 三种矩阵方程的解法：矩阵分解（相抵分解、满秩分解、奇异值分解）；Penrose定理求解；直积和矩阵拉直求解。

3. <font color='red'>**方法一：利用矩阵分解求解**</font>

   ==定理1.2.1==考察矩阵方程$AXB=D$，如果存在非奇异矩阵（**指行列式不为零的方阵，具有可逆性**）$P_i$和$Q_i$（$i=A,B$）满足式
   $$
   P_AAQ_A=\begin{bmatrix}
   I_r & 0 \\
   0 & 0
   \end{bmatrix},P_BBQ_B=\begin{bmatrix}
   I_r & 0 \\
   0 & 0
   \end{bmatrix}
   $$
   则方程$AXB=D$相容的充要条件是
   $$
   P_ADQ_B=\begin{bmatrix}
   J_{11} & 0 \\
   0 & 0
   \end{bmatrix}
   $$
   此时方程$AXB=D$的解为$X=Q_A\begin{bmatrix}
   J_r & Y_{12} \\
   Y_{21} & Y_{22}
   \end{bmatrix}P_B$。

   式中，$J_{11}\in\mathbb{C}^{r×r}$；$Y_{12},Y_{21},Y_{22}$是适当阶数的任意矩阵。由于$P_A、D、Q_B$是三个已知的矩阵，所以$J_{11}$可以计算出来。

   > [!NOTE]
   >
   > 证明：将$P_AAQ_A=\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix},P_BBQ_B=\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}$代入方程$AXB=D$可得
   > $$
   > P_A^{-1}\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}Q_A^{-1}XP_B^{-1}\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}Q_B^{-1}=D\\
   > \Rightarrow \\
   > \begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}Q_A^{-1}XP_B^{-1}\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}=P_ADQ_B
   > $$
   > 定义$Y=Q_A^{-1}XP_B^{-1}$，$J=P_ADQ_B$并对其分块
   > $$
   > Y=\begin{bmatrix}
   > 	Y_{11} & Y_{12} \\
   > Y_{21} & Y_{22}
   > \end{bmatrix},
   > J=\begin{bmatrix}
   > J_{11} & J_{12} \\
   > J_{21} & J_{22}
   > \end{bmatrix}
   > $$
   > 式中，$Y_{11}$和$J_{11}\in\mathbb{C}^{r×r}$。由此，式（12）简化为
   > $$
   > \begin{bmatrix}
   > Y_{11} & 0 \\
   > 0 & 0
   > \end{bmatrix}=\begin{bmatrix}
   > J_{11} & J_{12} \\
   > J_{21} & J_{22}
   > \end{bmatrix}=P_ADQ_B
   > $$
   > 所以，方程$AXB=D$相容当且仅当
   > $$
   > P_ADQ_B=\begin{bmatrix}
   > J_{11} & 0 \\
   > 0 & 0
   > \end{bmatrix}
   > $$
   > $P_ADQ_B$是三个已知矩阵的相乘，所以可以计算出$J_{11}$。
   >
   > 此时注意到$Y=Q_A^{-1}XP_B^{-1}$可得矩阵方程的解为
   > $$
   > X=Q_A\begin{bmatrix}
   > J_{11} & Y_{12} \\
   > Y_{21} & Y_{22}
   > \end{bmatrix}P_B
   > $$
   > 式中$J_{11}\in\mathbb{C}^{r×r}$，$Y_{12},Y_{21}$和$Y_{22}$是适当阶数的任意矩阵。

   > [!TIP]
   >
   > 上面的推导可以帮助理解，从解题的角度，我们只需要记住以下几个步骤：
   >
   > Step 1：$P_ADQ_B=\begin{bmatrix}
   > J_{11} & 0 \\
   > 0 & 0
   > \end{bmatrix}$（相容性检查）；
   >
   > Step 2：方程的解$X=Q_A\begin{bmatrix}
   > J_{11} & Y_{12} \\
   > Y_{21} & Y_{22}
   > \end{bmatrix}P_B$。

   > [!CAUTION]
   >
   > 若$r_A=r_B=r$，由定理1.2.1知，方程$AXB=D$相容的必要条件为$r_D\leq r$。也就是说，从$r_D\leq r$可以推导出方程$AXB=D$相容。（因为$D$是由三个矩阵相乘得来，他的秩不可能高于原矩阵的秩。相等的时候分解唯一。）

4. ==推论1.2.1==考察矩阵方程$AXA=A$，存在非奇异矩阵$P$和$Q$使得$PAQ=\begin{bmatrix}
   I_r & 0 \\
   0 & 0
   \end{bmatrix}$，则$G$是该方程解的充要条件是
   $$
   G=Q\begin{bmatrix}
   I_r & K \\
   L & M
   \end{bmatrix}P
   $$
   式中，$K,L,M$是适当阶数的任意矩阵，当$A$为可逆矩阵时，$G$唯一。

   > [!CAUTION]
   >
   > 推论1.2.1表明矩阵方程$AXA=A$的解一定存在$\Rightarrow$任一矩阵的减号逆一定存在。（可以出判断题）

5. ==例1.2.1==求解矩阵方程$AXA=A$，其中$A=\begin{bmatrix}
   1&0&1 \\
   0&1&1
   \end{bmatrix}$。`//TODO`

   > [!IMPORTANT]
   >
   > 这一题目的解题过程和解析如下：
   >
   > 解：$\begin{bmatrix}
   > A&I \\
   > I&*
   > \end{bmatrix}=\begin{bmatrix}
   > 1&0&1&1&0 \\
   > 0&1&1&0&1 \\ 1&0&0&*&* \\ 0&1&0&*&* \\ 0&0&1&*&*
   > \end{bmatrix}\Rightarrow(初等变换)\begin{bmatrix}
   > 1&0&0&1&0 \\
   > 0&1&0&0&1 \\ 1&0&-1&*&* \\ 0&1&-1&*&* \\ 0&0&1&*&*
   > \end{bmatrix}$，其中$*$为任意复数。
   >
   > 初等变换包括行变换和列变换，对于可逆的初等变换矩阵，左乘为行变换，右乘为列变换。
   >
   > 由此解得$P=\begin{bmatrix}
   > 1&0 \\
   > 0&1 
   > \end{bmatrix}$，$Q=\begin{bmatrix}
   > 1&0&-1 \\
   > 0&1&-1 \\ 0&0&1
   > \end{bmatrix}$。故方程$AXA=A$的解可以表示为
   > $$
   > G=Q\begin{bmatrix}
   > 1&0 \\
   > 0&1 \\ a&b
   > \end{bmatrix}P=\begin{bmatrix}
   > 1-a&-b \\
   > -a&1-b \\ a&b
   > \end{bmatrix}
   > $$
   > 式中$a,b$为任意数。

   > [!CAUTION]
   >
   > 推论1.2.1中可逆矩阵$P$和$Q$可通过初等变换求解，即
   > $$
   > \begin{bmatrix}
   > P&0 \\
   > 0&I 
   > \end{bmatrix}
   > \begin{bmatrix}
   > A&I \\
   > I&* 
   > \end{bmatrix}
   > \begin{bmatrix}
   > Q&0 \\
   > 0&I 
   > \end{bmatrix}=
   > \begin{bmatrix}
   > PAQ&P \\
   > Q&* 
   > \end{bmatrix}
   > $$
   > **为什么？**（初等变换矩阵均为可逆矩阵，而且根据推论1.2.1，$PAQ=\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}$，可以通过初等变换求解可逆矩阵$P$和$Q$）

6. ==定理1.2.2==考察矩阵方程$AXB=D$，若矩阵$A$和$B$有奇异值分解
   $$
   A=U_A\begin{bmatrix}
   \Sigma_r & 0 \\
   0 & 0
   \end{bmatrix}V_A^H,B=U_B\begin{bmatrix}
   \Lambda_r & 0 \\
   0 & 0
   \end{bmatrix}V_B^H
   $$
   则方程$AXB=D$相容的充分必要条件是
   $$
   U_A^HDV_B=\begin{bmatrix}
   J_{11} & 0 \\
   0 & 0
   \end{bmatrix}
   $$
   此时方程的解为
   $$
   X=V_A\begin{bmatrix}
   \Sigma_r^{-1}J_{11}\Lambda_r^{-1} & Y_{12} \\
   Y_{21} & Y_{22}
   \end{bmatrix}U_B^H
   $$

7. ==推论1.2.2==考察矩阵方程$AXA=A$，若矩阵$A$有奇异值分解$A=U\begin{bmatrix}
   \Sigma_r^{-1} & 0 \\
   0 & 0
   \end{bmatrix}V^H$，则矩阵$G$是该方程解的充要条件是
   $$
   G=V\begin{bmatrix}
   \Sigma_r^{-1} & K \\
   L & M
   \end{bmatrix}U^H
   $$
   式中$\Sigma_r$是$r$阶对角矩阵，其对角线元素是矩阵$A$的$r$个正奇异值；$K,L,M$是适当阶数的任意矩阵。当$A$为可逆矩阵时，$G$唯一。

8. <font color='red'>**方法二：利用Penrose定理求解**</font>

   ==定理1.2.3（Penrose定理）==矩阵方程$AXB=D$相容的充要条件为$AA^{(1)}DB^{(1)}B=D$，其中$A^{(1)}\in A\{1\}$，$B^{(1)}\in B\{1\}$；此时方程$AXB=D$的通解为
   $$
   X=A^{(1)}DB^{(1)}+Y-A^{(1)}AYBB^{(1)}
   $$
   式中，$Y$为任一$n×p$矩阵。
   $$
   AA^{(1)}A=A,BB^{(1)}B=B
   $$
   找到一个减号逆，遍历矩阵$Y$，即可构造出减号逆的集合$X=A^{(1)}DB^{(1)}+Y-A^{(1)}AYBB^{(1)}$。

   > [!NOTE]
   >
   > Penrose定理证明：
   >
   > 1、充分性（条件成立推结论）：设$AA^{(1)}DB^{(1)}B=D$成立，则参考矩阵方程$AXB=D$的形式，有$X=A^{(1)}DB^{(1)}$使得矩阵方程$AXB=D$相容，即充分性成立；
   >
   > 2、必要性（结论成立推条件）：假设矩阵方程$AXB=D$有解，不妨设$X$是其中一个解。根据推论1.2.1，任何一个矩阵都有其减号逆，则对于矩阵$A$、$B$来说，有$AA^{(1)}A=A$，$BB^{(1)}B=B$，代入原方程可得
   > $$
   > AA^{(1)}AXBB^{(1)}B=D
   > $$
   > 因为$AXB=D$，所以有$AA^{(1)}DB^{(1)}B=D$。必要性成立。
   >
   > 通解的证明：
   >
   > 1、充分性（条件成立推结论）：假设条件成立，即通解为定理中的形式，则将原公式左侧代入，得
   > $$
   > AXB=AA^{(1)}DB^{(1)}B+AYB-AA^{(1)}AYBB^{(1)}B=\\ AA^{(1)}DB^{(1)}B+AYB-AYB=AA^{(1)}DB^{(1)}B=D
   > $$
   > 通解形式的充分性成立。
   >
   > 2、必要性（结论成立推条件）：
   >
   > 假设$X$是矩阵方程$AXB=D$的任意解，则设
   > $$
   > X=A^{(1)}DB^{(1)}+X-A^{(1)}DB^{(1)}
   > $$
   > 将$AXB=D$代入上述公式，有
   > $$
   > X=A^{(1)}DB^{(1)}+X-A^{(1)}AXBB^{(1)}
   > $$
   > 符合通解形式。因此必要性成立。

9. ==推论1.2.3==非齐次线性方程组$Ax=b$相容的充要条件是$AA^{(1)}b=b$，其通解为
   $$
   x=A^{(1)}b+(I-A^{(1)}A)y, \forall y\in\mathbb{C}^n
   $$

10. ==推论1.2.4==齐次线性方程组$Ax=0$的通解为
    $$
    x=(I-A^{(1)}A)y, \forall y\in\mathbb{C}^n
    $$

11. ==推论1.2.5==设$A\in\mathbb{C}^{m×n}$，则
    $$
    A\{1\}=A^{(1)}+Z-A^{(1)}AZAA^{(1)}
    $$
    其中，$Z$为任一$n×m$矩阵。

    关于上述三个推论的例题：==例1.2.2==

12. <font color='red'>**方法三：利用直积和矩阵拉直求解**</font>

    ==定义1.2.2（Kronecker积）==设$A=(a_{ij})\in\mathbb{C}^{m×n}$，$B=(b_{ij})\in\mathbb{C}^{p×q}$，称如下<u>**分块矩阵**</u>
    $$
    A\otimes B=\begin{bmatrix}
    a_{11}B&a_{12}B&···&a_{1n}B \\
    a_{21}B&a_{22}B&···&a_{2n}B \\ \vdots&\vdots&&\vdots \\ a_{m1}B&a_{m2}B&···&a_{mn}B
    \end{bmatrix}\in\mathbb{C}^{mp×nq}
    $$
    为$A$与$B$的Kronecker积（或直积、张量积），简记为$A\otimes B=(a_{ij}B)$。

    辅助理解的例题：==例1.2.3==

13. ==命题1.2.1（直积的性质）==矩阵的Kronecker积具有以下性质：

    > [!IMPORTANT]
    >
    > （1）对任意复数 $k$，$(kA) \otimes B = A \otimes (kB) = k(A \otimes B)$;
    >
    > （2）$(A \otimes B) \otimes C = A \otimes (B \otimes C)$;
    >
    > （3）$A \otimes (B + C) = A \otimes B + A \otimes C$;
    >
    > （4）$(A \otimes B)^H = A^H \otimes B^H$;
    >
    > （5）若矩阵 $A$ 和 $C$，矩阵 $B$ 和 $D$ 均可相乘，则 $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$;
    >
    > （6）$\operatorname{rank}(A \otimes B) = \operatorname{rank}(A)\operatorname{rank}(B)$;
    >
    > （7）$(A \otimes B)^+ = A^+ \otimes B^+$.

14. ==定理1.2.4==设$f(x,y)=\sum_{i,j=0}^{k}c_{ij}x^iy^j$是变量$x,y$的复二元多项式，对任意矩阵$A\in\mathbb{C}^{m×m}$和$B\in\mathbb{C}^{n×n}$，定义矩阵
    $$
    f(A,B)=\sum_{i,j=0}^{k}c_{ij}(A^i\otimes  B^j)\in\mathbb{C}^{mn×mn}
    $$

    式中，$A^0=I_m$，$B^0=I_n$。若$A$和$B$的特征值分别为$\lambda_1,...,\lambda_m$和$\mu_1,...,\mu_n$，则$f(A,B)$的特征值为$f(\lambda_i,\mu_j)$，$i=1,...,m$，$j=1,...,n$。

15. ==推论1.2.6==设$A\in\mathbb{C}^{m×m}$，$B\in\mathbb{C}^{n×n}$，其特征值分别为$\lambda_1,...,\lambda_m$和$\mu_1,...,\mu_n$，则有

    （1）$A\otimes B$的特征值为$\lambda_i\mu_j$，$i=1,...,m$，$j=1,...,n$。

    （2）$A\otimes I_n+I_m\otimes B$的特征值为$\lambda_i+\mu_j$，$i=1,...,m$，$j=1,...,n$。

16. ==定义1.2.3（列拉直）==设$A=(a_{ij})\in\mathbb{C}^{m×n}$，并记$\boldsymbol{a_i}=[a_{1i},a_{2i},...,a_{mi}]^T$（$\boldsymbol{a_i}$此时是一个列向量），$i=1,2,...,n$。令
    $$
    vec(A)=\begin{bmatrix}
    \boldsymbol{a_1} \\
    \boldsymbol{a_2} \\
    \vdots \\
    \boldsymbol{a_n}
    \end{bmatrix}
    $$
    称为矩阵$A$的列拉直（或列展开）。（$vec(A)$是一个列向量，说白了就是将矩阵$A$的每一列前后接起来形成一个新的列向量。）

17. ==定理1.2.5==设$A\in\mathbb{C}^{m×n}$，$B\in\mathbb{C}^{n×p}$和$C\in\mathbb{C}^{p×q}$，则
    $$
    vec(ABC)=(C^T\otimes A)vec(B)
    $$

18. ==定理1.2.6==矩阵方程$AXB=D$相容的充要条件为
    $$
    rank(B^T\otimes A,vec(D))=rank(B^T\otimes A)
    $$
    定理1.2.3和定理1.2.6分别给出了矩阵方程$AXB=D$相容的充要条件。显然，这两个条件应该是互为等价的。

    这里可以理解为增广矩阵的秩等于系数矩阵的秩。这里是针对矩阵方程的推广。

19. ==应用==Lyapunov方程，定义为$AX+XA^H=-Q$，该方程在控制理论、通讯和动力系统中起着非常重要的作用.。我们常根据Lyapunov矩阵方程的解来检测系统的稳定性、可控性和可观测性等问题。

20. ==定理1.2.7==设$A \in \mathbb{C}^{n \times n}$，$Q \in \mathbb{C}^{n \times n}$，若矩阵$A$的所有特征值均具有负实部，则矩阵方程$AX + XA^H = -Q$ 有唯一解，且解$X \in \mathbb{C}^{n \times n}$ 可表示为
    $$
    X=\int_{0}^{+\infty}e^{At}Qe^{A^H t}dt
    $$

    进一步，若 \(Q\) 是（半正定、正定）Hermite矩阵，则解$X$也是（半正定、正定）Hermite矩阵.

### 1.1.3 减号逆

1. 推论1.2.1和推论1.2.2给出了矩阵减号逆的存在性结论，并指出任一矩阵的减号逆必存在。定理1.2.3也提供了矩阵减号逆的通解。此外，我们还可以利用**满秩分解**给出减号逆的一个特解。特别地，若$A\in\mathbb{C}_n^{m×n}$，则$A$的左逆$(A^HA)^{-1}A^H\in A\{1\}$；若$A\in\mathbb{C}_m^{m×n}$，则$A$的右逆$A^H(AA^H)^{-1}\in A\{1\}$。

   （矩阵列满秩有左逆，行满秩有右逆。）

   > [!NOTE]
   >
   > 为什么呢？因为矩阵列满秩的时候，有$m\geq n$，则$A^TA_{n}$为可逆矩阵，有$(A^TA)^{-1}A^TA=I$。所以$A$列满秩时，左逆为$(A^TA)^{-1}A^T$。行满秩右逆同理。

2. ==注1：==减号逆不唯一。若$A\in\mathbb{C}^{m×n}_r$，则其减号逆阶数为$n×m$（指矩阵的尺寸）。

3. ==例1.3.1==求矩阵$A=\begin{bmatrix}
   1&0&1 \\
   0&1&1
   \end{bmatrix}$的一个减号逆。

   > [!IMPORTANT]
   >
   > 可以看到矩阵$A$为行满秩矩阵，则矩阵有右逆$A^T(AA^T)^{-1}\in A\{1\}$，为$A$的一个减号逆。
   >
   > 为$A=\begin{bmatrix}
   > \frac{2}{3}&-\frac{1}{3} \\ 
   > -\frac{1}{3}&\frac{2}{3} \\ \frac{1}{3}&\frac{1}{3}
   > \end{bmatrix}$。计算二阶矩阵逆的方法：主对调负反号，除以行列式。

4. ==注2：==尽管利用奇异值分解或满秩分解求解矩阵的减号逆表达比较简洁，但利用相抵分解求解的计算量相对较小。

5. ==命题1.3.1（减号逆的性质）==

   > [!IMPORTANT]
   >
   > （1）$rank(A)\leq rank(A^-)$；
   >
   > （2）$AA^-$和$A^-A$均为幂等矩阵，且$rank(AA^-)=rank(A^-A)=rank(A)$；
   >
   > ​	$rank(A)\geq rank(AA^-)\geq rank(AA^-A)=rank(A)$；
   >
   > （3）若$B=P^{-1}AQ^{-1}$，则$QA^-P\in B\{1\}$。

6. ==【应用】==考察非齐次线性方程组$A\boldsymbol{x}=\boldsymbol{b}$。式中，$A\in\mathbb{C}^{m×n}$和$\boldsymbol{b}\in\mathbb{C}^m$给定，$\boldsymbol{x}\in\mathbb{C}^n$为待定向量。

   首先讨论$A^{(1)}$与方程组$A\boldsymbol{x}=\boldsymbol{b}$相容性的关系。

7. ==命题1.3.2（线性方程组相容条件）==设$A\in\mathbb{C}^{m×n}$，$\boldsymbol{b}\in\mathbb{C}^m$，考察非齐次线性方程组$A\boldsymbol{x}=\boldsymbol{b}$，则以下表达式等价：

   > [!IMPORTANT]
   >
   > （1）$rank(A)=rank(A,\boldsymbol{b})$
   >
   > （2）$\boldsymbol{b}\in R(A)$
   >
   > （3）$AA^{(1)}\boldsymbol{b}=\boldsymbol{b}$（Penrose定理）

8. ==定理1.3.1（减号逆与方程解的关系）==设$A\in\mathbb{C}^{m×n}$，则$G\in A\{1\}$的充要条件是$\boldsymbol{x}=G\boldsymbol{b}$是相容方程组$A\boldsymbol{x}=\boldsymbol{b}$的解。`//TODO` 找一下证明

9. ==推论1.3.1（减号逆与方程解的关系）==设$A\in \mathbb{C}^{m×n}$，非零向量$\boldsymbol{b}\in\mathbb{C}^m$，则相容方程组$A\boldsymbol{x}=\boldsymbol{b}$的解一定可以用$\boldsymbol{x}=G\boldsymbol{b}$表示，其中$G\in A\{1\}$。

10. ==例1.3.2==考察$\boldsymbol{x}=A^{(1)}\boldsymbol{b}$与相容方程组$A\boldsymbol{x}=\boldsymbol{b}$解的关系，其中$A=\begin{bmatrix}
    1&0&0 \\
    0&1&0
    \end{bmatrix}$，$\boldsymbol{b}=\begin{bmatrix}
    1 \\
    0
    \end{bmatrix}$。

    > [!IMPORTANT]
    >
    > Step（1）由于$\boldsymbol{b}\in R(A)$，故方程组$A\boldsymbol{x}=\boldsymbol{b}$相容。（命题1.3.2，线性方程组相容条件）
    >
    > Step（2）通过求解线性方程组（常规求法），可得该方程组的解为$\boldsymbol{x}=\begin{bmatrix}
    > 1 \\
    > 0 \\ c
    > \end{bmatrix}$，其中$\forall c\in \mathbb{C}$。
    >
    > Step（3）由减号逆的定义（直接把$\boldsymbol x$的四个元素设出来计算即可），可以求解得到$A^{(1)}=\begin{bmatrix}
    > 1&0 \\
    > 0&1 \\ a&b
    > \end{bmatrix}$，$\forall a,b\in \mathbb{C}$。则有$A^{(1)}\boldsymbol{b}=\begin{bmatrix}
    > 1 \\
    > 0 \\ a
    > \end{bmatrix}$，$\forall a\in\mathbb{C}$。可以看到，该解的形式和$\boldsymbol{x}$一模一样。那么$\boldsymbol{x}=A^{(1)}\boldsymbol b$。
    >
    > 另外，若$\boldsymbol b=0$，则不论$A^{(1)}$是任何矩阵，方程的解都是$0$。此时，齐次方程组的解为$\boldsymbol x=[0,0,c]^T,\forall c\in\mathbb{C}$。相容方程组的解不能用$\boldsymbol x=A^{(1)}\boldsymbol b$来表示。

11. ==定理1.3.2（相容线性方程组唯一解条件）==设$A\in\mathbb{C}^{m×n}$，$A^{(1)}A=I_n$当且仅当矩阵$A$是列满秩的。此时，相容方程组$A\boldsymbol{x}=\boldsymbol{b}$有唯一解。（其实就是以前学过的列满秩系数矩阵构成的线性方程组具有唯一解）

    > [!NOTE]
    >
    > 证明：此处证的是当且仅当的部分。
    >
    > 充分性（条件成立推结论）：由命题1.3.1减号逆的性质第二条知，$A^{(1)}A$是幂等矩阵。根据幂等矩阵的性质，它的特征值只能是$0$或$1$。又因为$A$是列满秩矩阵，所以特征值只能全部是$1$。这样，$A^{(1)}A$相似于单位矩阵，即$A^{(1)}A=I_n$。
    >
    > 必要性（结论成立推条件）：设$A^{(1)}A=I_n$，矩阵$A$有左逆$A^{(1)}$，根据**定理3.2.3（教材P99定理+证明）**，矩阵$A$有左逆的充要条件是矩阵$A$是列满秩矩阵。

12. ==注3：==当矩阵$A$是列满秩矩阵时，它的任一减号逆$A^{(1)}$都是$A$的左逆。当矩阵$A$是行满秩矩阵时，它的任一减号逆$A^{(1)}$都是$A$的右逆。

13. ==例1.3.3==求向量$\boldsymbol b\in\mathbb{C}^m$在$R(A)$上的正交投影，其中$A\in\mathbb{C}^{m×n}$。（教材P248 例5.3.3）

    这里几个概念要先捋清楚。

    ==正交投影矩阵：==回顾定义3.7.2（幂等矩阵），幂等矩阵也称为投影矩阵。Hermite幂等矩阵称为正交投影矩阵。

    在教材P131的例3.7.2，矩阵$A$是列满秩矩阵，可以设$P=A(A^HA)^{-1}A^H$，$P^2=P$，为幂等矩阵。

    但是在本例中，并没有说明矩阵$A$是不是列满秩矩阵，因此$(A^HA)^{-1}$不一定存在。

    我们把$(A^HA)^{-1}$更换为减号逆$(A^HA)^-$，并令$P=A^H(A^HA)^-A$。

14. ==定理1.3.3==对任意矩阵$A\in\mathbb{C}^{m×n}$，$A(A^HA)^{-}A^H$是正交投影矩阵。

15. ==注4：==对一般的矩阵$A$，$(AA^H)^-=(A^H)^-A^-$。

16. ==推论1.3.2==对任意矩阵$A\in\mathbb{C}^{m×n}$，$A(A^HA)^-A^HA=A$。（证明见教材P249 推论5.3.2）。

    > [!TIP]
    >
    > 另一种证明方法：
    >
    > 设 $G = (A^H A)^-$，已知 $A^H A G A^H A = A^H A$。
    >
    > 考虑 $A - A G A^H A$，左乘 $A^H$：
    >
    > $$
    > A^H (A - A G A^H A) = A^H A - A^H A G A^H A = 0
    > $$
    >
    > 所以 $R(A - A G A^H A) \subseteq N(A^H)$。
    >
    > 但显然又有 $R(A - A G A^H A) \subseteq R(A)$。
    >
    > 而 $R(A) \cap N(A^H) = \{0\}$（这是基本空间的正交补关系），因此：
    >
    > $$
    > A - A G A^H A = 0 \implies A(A^H A)^- A^H A = A
    > $$

17. ==例1.3.4（最小二乘问题，教材P249 例5.4.4（书上序号打错了，应该是5.3.4））==考察线性方程组$A\boldsymbol x=\boldsymbol b$，其中，矩阵$A\in\mathbb{C}^{m×n}$和向量$\boldsymbol b\in\mathbb{C}^m$给定，向量$\boldsymbol x\in\mathbb{C}^n$待定。求向量$\boldsymbol x$使得$||A\boldsymbol x-\boldsymbol b||_2$最小。

    最小二乘问题是指$b$在$A$的列空间上的正交投影。

### 1.1.4 极小范数广义逆

1. ==定理1.4.1==（教材P250 定理5.4.1，无证明）设矩阵$A\in\mathbb{C}^{m×n}$的奇异值分解为
   $$
   A=U\begin{bmatrix}
   \Sigma_r & 0 \\
   0 & 0
   \end{bmatrix}V^H
   $$
   则$G\in A\{1,4\}$的充要条件是
   $$
   G=V\begin{bmatrix}
   \Sigma_r^{-1} & K \\
   0 & M
   \end{bmatrix}U^H
   $$
   式中，$K$和$M$为适当阶数的任意矩阵。

2. ==注1：==若$A\in\mathbb{C}^{m×n}_n$，则$(A^HA)^{-1}A^H\in A\{1,4\}$；若$A\in\mathbb{C}^{m×n}_m$，则$A^H(AA^H)^{-1}\in A\{1,4\}$。

3. 【思考】能否利用Penrose定理求解矩阵的极小范数广义逆？第一个条件可以使用Penrose定理，但是第四个条件$(XA)^H=XA$很难使用Penrose定理求解。为了解决该问题，引入定理1.4.2。

4. ==定理1.4.2==（教材P250 定理5.4.2，证明仍然从充分性和必要性两个方向来证明）设$A\in\mathbb{C}^{m×n}$，则
   $$
   A\{1,4\}=\{X\in\mathbb{C}^{n×m}|XA=A_m^-A\}
   $$
   该定理是指，$A\{1,4\}$满足方程$XA=A\{1,4\}A$。

5. ==定理1.4.3==（教材P250 定理5.4.3）设$A_m^-\in A\{1,4\}$，则
   $$
   A\{1,4\}=\{A_m^-+Z(I-AA_m^-)|Z\in\mathbb{C}^{n×m}\}
   $$

   > [!IMPORTANT]
   >
   > 证明：
   >
   > 由Penrose定理可得，方程$XA=A_m^-A$的通解为
   > $$
   > X=A_m^-AA_m^-+Y-YAA_m^-,\forall Y\in\mathbb{C}^{n×m}
   > $$
   > 设$Z$为$\mathbb{C}^{n×m}$中任意矩阵，将$Y=A_m^-+Z$代入上式可得
   > $$
   > X=A_m^-+Z(I-AA_m^-)
   > $$

6. 【应用】当相容矩阵$A\boldsymbol x=\boldsymbol b$有多个解时，往往需要从中挑选一个合适的解。此时我们一般挑选所有解中范数最小的解，或者是其中之一，即如下优化问题
   $$
   \mathop{\mathop{min}\limits_{A\boldsymbol x=\boldsymbol b}}\limits_{AA^{(1)}\boldsymbol b=\boldsymbol b}||\boldsymbol x||_2
   $$
   将满足该优化问题的解称为相容线性方程组$A\boldsymbol x=\boldsymbol b$的极小范数解。

7. ==定理1.4.4==（）





### 1.1.5 最小二乘广义逆





### 1.1.6 加号逆





## 1.2 矩阵微积分及其应用







## 1.3 投影矩阵及其应用









# 2.题目与解法

