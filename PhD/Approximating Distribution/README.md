# Approximating the Invariant Wealth Distribution
## 1. Preliminaries
This is an exercise to implement a list of different methods to approximate the invariant wealth distribution of ABH (Aiyagari-Bewley-Huggett) type models. The methods to the implemented are the following:
- Eigen Vector Method
- Monte Carlo Simulation
- Discretization of the CDF
- Discretization of the PDF
- Piece Wise Linear interpolation (Still to be done) 
- Collocation (Still to be done)

The methodology of the listed approximation methods are described in detail in:
- [Lecture Notes Gianluca Violante: Aproximating the Invariant distribution (2015)](http://www.econ.nyu.edu/user/violante/NYUTeaching/QM/Fall15/Lectures/Lecture11_Distribution_Slides.pdf). 
- [QuantEcon Website](https://lectures.quantecon.org/py/aiyagari.html) of Thomas,Sargent and  John,Stachurski
- [A method for solving and estimating heterogeneous agent macro models, Winberry (2018)](https://onlinelibrary.wiley.com/doi/pdf/10.3982/QE740)
## 2.The Model
### Households
Infinitely lived households solve the following maximization problem:

<img src="https://render.githubusercontent.com/render/math?math=\max \Big\{ \mathbb{E} \sum_{t=0}^{\infty} \beta^{t}u(c_{t}) \Big\}">

Subject to: 

<img src="https://render.githubusercontent.com/render/math?math=a_{t\pm1}\pm c_{t}\leq wz_{t}\pm (1\pmr)a_{t}">

<img src="https://render.githubusercontent.com/render/math?math=c_{t}\geq0 \,\,\,\,\mbox{and}\,\,\,\,a_{t}\geq0">
Where:

- <img src="https://render.githubusercontent.com/render/math?math=c_{t}"> is current consumption.
- <img src="https://render.githubusercontent.com/render/math?math=a_{t}"> is assets.
- <img src="https://render.githubusercontent.com/render/math?math=z_{t}"> exogenous component of labor income.
- <img src="https://render.githubusercontent.com/render/math?math=w"> is the wage rate.
- <img src="https://render.githubusercontent.com/render/math?math=r"> is the interest rate.
- Agents are not allowed to borrow.

The Exogenous process <img src="https://render.githubusercontent.com/render/math?math=z_{t}"> follows a finite Markov Chain with given stochastic Matrix. The wage rate and interest rate are given. Households supply labor inelastically.

### Firms
Firms produce output by hiring capital and labor, with constant returns to scale technology:

<img src="https://render.githubusercontent.com/render/math?math=Y_{t}=AK_{t}^{\alpha}N^{1-\alpha}">

With:

- <img src="https://render.githubusercontent.com/render/math?math=A=1"> Total Factor Productivity.
- <img src="https://render.githubusercontent.com/render/math?math=\alpha"> Capital share.
- <img src="https://render.githubusercontent.com/render/math?math=K_{t}"> is Aggregate Capital.
- <img src="https://render.githubusercontent.com/render/math?math=N"> is Total labor supply.

Then the firms maximization problem is:

<img src="https://render.githubusercontent.com/render/math?math=\max_{K,N} \Big[AK_{t}^{\alpha}N^{1-\alpha} -(r\pm\delta)K-w N\Big]">

Where <img src="https://render.githubusercontent.com/render/math?math=\delta"> is the depreciation rate

The code is written in:
```
Julia 1.0.2
```
