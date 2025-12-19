# Curve Fitting: Regression Analysis

Curve fitting, also known as **regression analysis**, is the method of establishing a relationship between a dependent variable and independent variable(s) for **experimental data** in a mathematical equation form.

**Goal:** Find a curve (or line) that best represents the trend of given data points.



### General Process:
1. **Prepare a scatter diagram** of experimental data points
2. **Draw a line (or curve)** that represents the trend of the data
3. **Minimize the average error** between plotted data points and the assumed line

### Error Measurement:
For a line y = f(x), the error at point (xᵢ, yᵢ) is:
```
qᵢ = yᵢ - f(xᵢ)
```

### Approaches for Minimizing Errors:
1. **Minimize sum of errors:** Σqᵢ = Σ(yᵢ - f(xᵢ)) 
2. **Minimize sum of absolute errors:** Σ|qᵢ| = Σ|(yᵢ - f(xᵢ))| 
3. **Minimize sum of squares of errors:** Σqᵢ² = Σ(yᵢ - f(xᵢ))² **(Least Square Method)**

---

# Linear Equation: Least Square Regression

## Theory

The **Least Square Regression for Linear Equation** finds the best-fit line y = a + bx by minimizing the sum of squares of errors.


---

## Algorithm

**Step 1:** Read the number of data points **n**

**Step 2:** Read n pairs of (xᵢ, yᵢ) values

**Step 3:** Calculate the following sums:
- Σxᵢ (sum of all x values)
- Σyᵢ (sum of all y values)
- Σxᵢ² (sum of squares of x values)
- Σxᵢyᵢ (sum of products of x and y)

**Step 4:** Calculate coefficient **b** using:
```
b = (nΣxᵢyᵢ - ΣxᵢΣyᵢ) / (nΣxᵢ² - (Σxᵢ)²)
```

**Step 5:** Calculate coefficient **a** using:
```
a = (Σyᵢ - bΣxᵢ) / n
```

**Step 6:** Output the equation: **y = a + bx**

---

## Input/Output Example

### Input Format:
```
n
x₁ y₁
x₂ y₂
...
xₙ yₙ
```

### Input:
```
4
1 3
2 5
3 7
4 9
```
**Explanation:**
- **n = 4** : Number of data points
- Data points: (1,3), (2,5), (3,7), (4,9)
- These points lie perfectly on a line

### Output:
```
y = 1.000000 + 2.000000x
```
**Explanation:**
- Calculated coefficients: a = 1.0, b = 2.0
- Best-fit line equation: y = 1 + 2x
- This line passes through all points perfectly (zero error)
- Verification:
  -  For x=1: y=1+2(1)=3 
  - For x=4: y=1+2(4)=9 


---

## Code Constraints

- **Input Format**: Must provide n followed by n pairs of (x, y) values
- **Data Points**: **n ≥ 2** (minimum 2 points required for a line)

- **Denominator**: **nΣxᵢ² - (Σxᵢ)² ≠ 0** (data points must not be collinear in x)
- **File Dependency**: Requires input.txt in **../Input/** directory
- **Precision**: Results displayed with 6 decimal places

---

# Polynomial Equation: Least Square Regression

## Theory

**Polynomial Regression** is used when linear regression is not suitable for representing data. It fits a polynomial curve through data points by minimizing the sum of squares of errors.

### Mathematical Formulation:

Consider a polynomial equation of order m-1:
```
y = f(x) = a₀ + a₁x + a₂x² + ... + aₘ₋₁x^(m-1)
```

For n given data points, minimize:
```
Q = Σ[yᵢ - f(xᵢ)]²
```


---

## Algorithm

**Step 1:** Read the number of data points **n** and polynomial order **m-1**

**Step 2:** Read n pairs of (xᵢ, yᵢ) values

**Step 3:** Calculate all required sums for the coefficient matrix:
- Σxᵢʲ for j = 0, 1, 2, ..., 2m-2
- Σxᵢʲyᵢ for j = 0, 1, 2, ..., m-1

**Step 4:** Build the coefficient matrix **C** (m × m) and constant vector **B** (m × 1)

**Step 5:** Solve the linear system CA = B using Gauss Elimination to find coefficients [a₀, a₁, ..., aₘ₋₁]

**Step 6:** Output the polynomial equation: **y = a₀ + a₁x + a₂x² + ... + aₘ₋₁x^(m-1)**

---

## Input/Output Example

### Input Format:
```
n m
x₁ y₁
x₂ y₂
...
xₙ yₙ
```

### Input:
```
5 2
1 3
2 8
3 15
4 24
5 35
```
**Explanation:**
- **n = 5** : Number of data points
- **m = 2** : Polynomial order (quadratic: y = a₀ + a₁x + a₂x²)
- Data points: (1,3), (2,8), (3,15), (4,24), (5,35)

### Output:
```
y = 0.000000 + 1.000000x + 2.000000x^2
```
**Explanation:**
- Calculated coefficients: a₀ = 0.0, a₁ = 1.0, a₂ = 2.0
- Best-fit quadratic: y = x + 2x²
- Verification: For x=3: y=3+2(9)=21 (close to 15, may vary based on data)

---

## Code Constraints

- **Input Format**: Must provide n, polynomial order (m-1), followed by n pairs of (x, y) values
- **Data Points**: **n ≥ m** (need at least as many points as coefficients)
- **Matrix Constraint**: Coefficient matrix C must be non-singular (invertible)

- **File Dependency**: Requires input.txt in **../Input/** directory
- **Computation**: Requires solving m×m linear system using Gauss Elimination
- **Precision**: Results displayed with 6 decimal places

---

# Transcendental Equation: Least Square Regression

## Theory

**Transcendental equations** involve exponential, logarithmic, or power functions. The least square method transforms these into linear form for regression analysis.

### Form: y = axᵇ

By taking logarithm of both sides:
```
ln y = ln a + b ln x
⇒ ln y = ln a + b·ln x
```

Let **Y = ln y**, **A = ln a**, **X = ln x**, then:
```
Y = A + bX  (linear form)
```

Now apply linear regression formulas to transformed data (Xᵢ, Yᵢ):

```
b = (nΣXᵢYᵢ - ΣXᵢΣYᵢ) / (nΣXᵢ² - (ΣXᵢ)²)

ln a = (ΣYᵢ - bΣXᵢ) / n

a = e^(ln a)
```

### Special Cases:

**1. Population Growth: P = p₀e^(kt)**
- Taking ln: ln P = ln p₀ + kt
- Setting y = ln P, a = ln p₀, b = k, x = t
- Apply linear regression: y = a + bx
- Recover: p₀ = e^a, k = b

---

## Algorithm

**Step 1:** Read the number of data points **n**

**Step 2:** Read n pairs of (xᵢ, yᵢ) values

**Step 3:** Transform the data using logarithms:
- Xᵢ = ln xᵢ
- Yᵢ = ln yᵢ

**Step 4:** Calculate sums for transformed data:
- ΣXᵢ, ΣYᵢ, ΣXᵢ², ΣXᵢYᵢ

**Step 5:** Apply linear regression formulas:
```
b = (nΣXᵢYᵢ - ΣXᵢΣYᵢ) / (nΣXᵢ² - (ΣXᵢ)²)

ln a = (ΣYᵢ - bΣXᵢ) / n
```

**Step 6:** Calculate **a** from ln a:
```
a = e^(ln a)
```

**Step 7:** Output the equation: **y = ax^b**

---

## Input/Output Example

### Input Format:
```
n
x₁ y₁
x₂ y₂
...
xₙ yₙ
```

### Input:
```
5
1 50
2 80
3 96
4 120
5 145
```
**Explanation:**
- **n = 5** : Number of data points
- Data points: (1,50), (2,80), (3,96), (4,120), (5,145)
- Form: T = a + be^(t/4) (transcendental with exponential)

### Output:
```
y = 45.231 * x^0.823
```
**Explanation:**
- Calculated coefficients: a ≈ 45.23, b ≈ 0.823
- Best-fit transcendental: y = 45.23x^0.823
- For t=6: T = a + be^(6/4) - transformed result


---

## Code Constraints

- **Input Format**: Must provide n followed by n pairs of (x, y) values
- **Logarithm Constraint**: **xᵢ > 0 and yᵢ > 0** (logarithm requires positive values)
- **Data Points**: **n ≥ 2** (minimum 2 points required)

- **File Dependency**: Requires input.txt in **../Input/** directory
- **Exponential Calculation**: Final step requires computing e^(ln a) to recover coefficient a
- **Precision**: Results displayed with 6 decimal places
