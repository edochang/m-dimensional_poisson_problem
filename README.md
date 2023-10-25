# 3D and m-dimensions Poisson Problem

In section 7.1.1 of the textbook (https://www.cs.utexas.edu/users/flame/laff/alaff/chapter07-why-sparse-2d.html) we learned about a simple Partial Differential Equation (PDE), Poison's Equation.  What would this problem look like in 3-dimensions?

## 3-Dimensional Poisson Problem Derivation ##

![3D Poisson Problem](/markdown_assets/poisson_problem_1.png)

Let $\phi_i$ be the value of $f(x,y)$ at the mesh point $i$.  One can approximate

$\begin{array}{l}
x-\mathrm{dimension}:\;\frac{\partial^2 u\left(x,y,z\right)}{\partial x^2 }\approx \frac{u\left(x-h,y,z\right)-2u\left(x,y,z\right)+u\left(x+h,y,z\right)}{h^2 }\\
y-\mathrm{dimension}:\;\frac{\partial^2 u\left(x,y,z\right)}{\partial y^2 }\approx \frac{u\left(x,y-h,z\right)-2u\left(x,y,z\right)+u\left(x,y+h,z\right)}{h^2 }\\
z-\mathrm{dimension}:\;\frac{\partial^2 u\left(x,y,z\right)}{\partial z^2 }\approx \frac{u\left(x,y,z-h\right)-2u\left(x,y,z\right)+u\left(x,y,z+h\right)}{h^2 }
\end{array}$

so that

$-\frac{\partial^2 u}{{\partial x}^2 }-\frac{\partial^2 u}{{\partial y}^2 }-\frac{\partial^2 u}{{\partial z}^2 }=f\left(x,y,z\right)$

becomes

$\begin{array}{l}
\frac{-u\left(x-h,y,z\right)+2u\left(x,y,z\right)-u\left(x+h,y,z\right)}{h^2 }+\frac{-u\left(x,y-h,z\right)+2u\left(x,y,z\right)-u\left(x,y+h,z\right)}{h^2 }+\\
\frac{-u\left(x,y,z-h\right)+2u\left(x,y,z\right)-u\left(x,y,z+h\right)}{h^2 }=f\left(x,y,z\right)
\end{array}$

or, equivalently,

$\frac{-u\left(x-h,y,z\right)-u\left(x,y-h,z\right)-u\left(x,y,z-h\right)+6u\left(x,y,z\right)-u\left(x+h,y,z\right)-u\left(x,y+h,z\right)-u\left(x,y,z+h\right)}{h^2 }=f\left(x,y,z\right)$

If $\left(x,y,z\right)$ corresponds to the point  in a mesh where the interior points form an N x N x N grid cube, this translates to the system of linear equations

$-\upsilon_{i-N^2 } -\upsilon_{i-N} -\upsilon_{i-1} +6\upsilon_i -\upsilon_{i+1} -\upsilon_{i+N} -\upsilon_{i+N^2 } =h^2 \phi_i$

This can be rewritten as

$\upsilon_i =\frac{h^2 \phi +\upsilon_{i-N^2 } +\upsilon_{i-N} +\upsilon_{i-1} +\upsilon_{i+1} +\upsilon_{i+N} +\upsilon_{i+N^2 } }{6}$

### Rephrase as Solving $\mathit{\mathbf{A}}\upsilon ={\mathit{\mathbf{h}}}^2 \phi$ ###

All this insights can be put together into a systems of linear equations where $\phi_i =f\left(\chi_i ,\psi_j ,\zeta_k \right)\;$ if $\left(\chi_i ,\psi_j ,\zeta_k \right)$ is the point assoicated with the value $\upsilon_i$.  In matrix notation this becomes: $A\upsilon =h^2 \phi$.  The figure below visually describes the matrix notation for a mesh (3-dimensional) cube of 3 x 3 x 3, with 27 interior points.

![3D Poisson Problem - 3 x 3 x 3](/markdown_assets/poisson_problem_2.png)

### Describe The Sparsity Pattern of the Matrix A for 3-Dimensions ###

Comparing the above figure in the previous section to a mesh (2-dimensional - 3 x 3) grid:

![3D Poisson Problem - 2D Mesh Grid](/markdown_assets/poisson_problem_3.png)

We can see the $A$ matrix has increase in size ( $N^3$ x $N^3$ ) to account for the N x N x N (or $N^3$) interior mesh points of a mesh cube.  We also see an extra band of -1 coefficients in the upper and lower triangular parts of the matrix, which represent the extra dimension in the negative and positive directions of the 7-point stencil: $-\upsilon_{i-N^2 }$ and $-\upsilon_{i+N^2 }$.  We can also see more sparsity when the interior mesh points' borders a boundary, so one-two points of the 7-point stencil sits on the boundary.  See the figure below for a visual example and pattern.

![3D Poisson Problem - Sparsity Pattern Highlights](/markdown_assets/poisson_problem_4.png)

### Educated Guess for discretizing Poisson Problem in an arbitrary, m ###

The A matrix will have m numbers of sub and super diagonals (below and above the main diagonal).  The main diagonal will be amplified and take on the value of $2*m$ as per working out the math in the section above.  Additional sparsity will be created based on interior points impacted by edges / boundaries.  If I have 4-dimension stencil, I will observe more sparsity spread across the 1-dimension diagonals, then less but longer bands of sparsity in the 2-dimension diagonals, and again get less, but longer bands of sparsity in the 3-dimension diagonals.  The 4-dimension diagonals will have no sparse bands in its diagonals, because the boundary impacts exist outside of matrix A.  Looking closer at the length of the gaps in each dimension we can see it follows a pattern of $N^{\mathrm{power}}$.  1-dimension diagonal gap length is $N^0 =1$, 2-dimension diagonal gap length is $N^2$... and this pattern continues up to the $m-1$ dimension.  We can identify the intersection of ( $i-N^{\mathrm{power}}$ and $i+N^{\mathrm{power}}$ ) per dimension by doing remainder division ( $\mathrm{mod}\left(i,N^{\mathrm{power}} \;\right)$ ).  

See my Create_mdiml_Poisson_problem_nzA() implement for more information.
- Creates the sparse format of the matrix corresponding to the left hand side matrix, A, in the discretized Poisson problem, , in dimension with 
- Work for any positive integer, , not just  = 1, 2, and 3.

### Large A Reviews ###
Used the Create_mdiml_Poisson_problem_nzA() code to try different variation of dimensions to validate my guess on the sparsity pattern.

**N = 2; m = 5**
![3D Poisson Problem - N = 2; m = 5](/markdown_assets/poisson_problem_5.png)

**N = 2; m = 6**
![3D Poisson Problem - N = 2; m = 6](/markdown_assets/poisson_problem_6.png)

**N = 3; m = 3**
![3D Poisson Problem - N = 3; m = 3](/markdown_assets/poisson_problem_7.png)

**N = 3; m = 4**
![3D Poisson Problem - N = 3; m = 4](/markdown_assets/poisson_problem_8.png)

**N = 4; m = 3**
![3D Poisson Problem - N = 4; m = 3](/markdown_assets/poisson_problem_9.png)

**N = 4; m = 4**
![3D Poisson Problem - N = 4; m = 4](/markdown_assets/poisson_problem_10.png)