# Runge-Kutta Method (RK4)

## Theory

### Ordinary Differential Equations (ODEs)

An **Ordinary Differential Equation (ODE)** is an equation that contains one independent variable and one or more of its derivatives with respect to that variable. In mathematical terms, an ODE represents a relationship between an independent variable x, a dependent variable y, and its derivatives y', y'', ..., yn.



### Order of ODE:
The **order** of an ODE is defined as the order of the highest derivative that occurs in the equation.

**Example:** dy/dx = f(x, y) is a **first-order ODE**

### Initial Value Problem:
An ODE with an initial condition is called an **Initial Value Problem (IVP)**:
```
dy/dx = f(x, y)  with  y(x₀) = y₀
```

Where:
- **x₀** = initial value of independent variable
- **y₀** = initial value of dependent variable y(x₀)

### Runge-Kutta Method:

The **Runge-Kutta (RK) method** is a family of iterative methods used to approximate the solutions of ODEs. The **4th-order Runge-Kutta method (RK4)** is the most commonly used and provides an excellent balance between accuracy and computational efficiency.


### RK4 Formula:

To compute the next value yₙ₊₁ from the current value yₙ:

```
yₙ₊₁ = yₙ + (1/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

Where the intermediate slopes are:

```
k₁ = h·f(xₙ, yₙ)

k₂ = h·f(xₙ + h/2, yₙ + k₁/2)

k₃ = h·f(xₙ + h/2, yₙ + k₂/2)

k₄ = h·f(xₙ + h, yₙ + k₃)
```

**Parameters:**
- **h** = step size (increment in x)
- **f(x, y)** = the derivative function dy/dx
- **k₁** = slope at the beginning of interval
- **k₂** = slope at midpoint using k₁
- **k₃** = slope at midpoint using k₂ (more accurate)
- **k₄** = slope at end of interval using k₃


---

## Algorithm

**Step 1:** Start with initial conditions
- x₀ = initial value of independent variable
- y₀ = initial value of y(x₀)
- Target x value where we want to find y(x)

**Step 2:** Choose a step size **h**
- The step size h defines how much to increment x in each step
- Smaller h gives more accuracy but requires more computations
- Calculate number of steps: n = (x - x₀) / h

**Step 3:** For each step i = 0 to n-1, compute four intermediate slopes:

Calculate **k₁** (slope at beginning):
```
k₁ = h·f(xₙ, yₙ)
```

Calculate **k₂** (slope at midpoint using k₁):
```
k₂ = h·f(xₙ + h/2, yₙ + k₁/2)
```

Calculate **k₃** (slope at midpoint using k₂):
```
k₃ = h·f(xₙ + h/2, yₙ + k₂/2)
```

Calculate **k₄** (slope at end using k₃):
```
k₄ = h·f(xₙ + h, yₙ + k₃)
```

**Step 4:** Update y using weighted average of slopes:
```
yₙ₊₁ = yₙ + (k₁ + 2k₂ + 2k₃ + k₄) / 6
```

**Step 5:** Increment x:
```
xₙ₊₁ = xₙ + h
```

**Step 6:** Repeat Steps 3-5 until reaching target x value

**Step 7:** Output the final y value

---

## Input/Output Example

### Input Format:
```
x₀ y₀ x h
```

### Input:
```
0.0 1.0 2.0 0.2
```
**Explanation:**
- **x₀ = 0.0** : Initial value of x
- **y₀ = 1.0** : Initial value of y (i.e., y(0) = 1)
- **x = 2.0** : Target x value where we want to find y
- **h = 0.2** : Step size

### Output:
```
=== Runge Kutta Method ===
Initial Condition (x0, y0): (0.000000, 1.000000)
Target x: 2.000000
Step size h: 0.200000

Result y(2.000000) = 1.103639
```
**Explanation:**
- ODE being solved: dy/dx = (x - y)/2 with y(0) = 1
- Number of steps: n = (2.0 - 0.0) / 0.2 = 10 steps
- At each step, compute k₁, k₂, k₃, k₄ and update y
- Final result: y(2.0) ≈ 1.103639

---

## Code Constraints

- **Input Format**: Must provide exactly 4 values: x₀, y₀, x, h
- **Step Size**: **h > 0** (positive step size required)
- **Step Size Limit**: **h ≤ (x - x₀)** (step size cannot exceed interval)
- **Target Value**: **x > x₀** (target must be greater than initial value)
- **Number of Steps**: **(x - x₀) / h** should be a positive integer 

- **Function**: Currently hardcoded as f(x,y) = (x-y)/2 in code
- **File Dependency**: Requires input.txt in **../Input/** directory
- **Precision**: Results displayed with 6 decimal places

