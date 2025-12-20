# Numerical Differentiation using Newton's Interpolation Method

## Theory

**Numerical differentiation** is the process of finding the derivative of a function using numerical methods when analytical differentiation is difficult or impossible. Newton's interpolation method provides a way to approximate derivatives using a finite difference table constructed from discrete data points.


### Newton's Forward Difference Formula:

For equally spaced points with spacing h, define:
- **u = (x - x₀) / h** (normalized distance from first point)
- **Δʸ** = forward differences of y

The interpolating polynomial is:
```
y(x) = y₀ + uΔy₀ + [u(u-1)/2!]Δ²y₀ + [u(u-1)(u-2)/3!]Δ³y₀ + ...
```

### First Derivative Formula:

Differentiating with respect to x:
```
dy/dx = (1/h)[Δy₀ + ((2u-1)/2)Δ²y₀ + ((3u²-6u+2)/6)Δ³y₀ + ((4u³-18u²+22u-6)/24)Δ⁴y₀ + ...]
```

Where:
- **h** = interval spacing = (b - a) / n
- **u** = (p - a) / h (at point p where derivative is needed)
- **Δʸy₀** = y-th order forward difference at x₀

### Second Derivative Formula:

Differentiating twice:
```
d²y/dx² = (1/h²)[Δ²y₀ + (u-1)Δ³y₀ + ((11u²-23u+10)/12)Δ⁴y₀ + ...]
```

### Forward Difference Table:

The forward difference table is constructed as:
```
Δy₀ = y₁ - y₀
Δ²y₀ = Δy₁ - Δy₀ = y₂ - 2y₁ + y₀
Δ³y₀ = Δ²y₁ - Δ²y₀
...
```

### Key Concepts:

1. **Forward Interpolation**: Used when the point p is near the beginning of the data range
2. **Backward Interpolation**: Used when the point p is near the end of the data range
3. **Difference Table**: Systematic arrangement of successive differences for computation

---

## Algorithm

**Step 1:** Read input parameters:
- **b** = upper limit
- **a** = lower limit
- **n** = number of data points
- **p** = point where derivative is to be calculated

**Step 2:** Calculate the step size:
```
h = (b - a) / n
```

**Step 3:** Generate data points:
- For i = 0 to n-1:
  - Xᵢ = a + i·h
  - Yᵢ = f(Xᵢ) (evaluate function at Xᵢ)

**Step 4:** Construct the forward difference table:
- Initialize first column: table[i][0] = Yᵢ
- For each order j = 1 to n-1:
  - For each row i = 0 to n-j-1:
    - table[i][j] = table[i+1][j-1] - table[i][j-1]

**Step 5:** Calculate normalized parameter:
```
u = (p - a) / h
```

**Step 6:** Extract differences from table:
```
Δy₀ = table[0][1]
Δ²y₀ = table[0][2]
Δ³y₀ = table[0][3]
Δ⁴y₀ = table[0][4]
```

**Step 7:** Calculate first derivative using forward formula:
```
dy/dx = [Δy₀ + ((2u-1)/2)Δ²y₀ + ((3u²-6u+2)/6)Δ³y₀ + ((4u³-18u²+22u-6)/24)Δ⁴y₀] / h
```

**Step 8:** Calculate second derivative:
```
d²y/dx² = [Δ²y₀ + (u-1)Δ³y₀ + ((11u²-23u+10)/12)Δ⁴y₀] / h²
```

**Step 9:** Calculate errors (if true values are known):
```
Error₁ = |((true value - calculated) / true value)| × 100%
Error₂ = |((true value - calculated) / true value)| × 100%
```

**Step 10:** Output results:
- Difference table
- First derivative (calculated and true value)
- Second derivative (calculated and true value)
- Percentage errors

---

## Input/Output Example

### Input Format:
```
b a n p
```

### Input:
```
1.0 0.0 11 0.1
```
**Explanation:**
- **b = 1.0** : Upper limit of interval
- **a = 0.0** : Lower limit of interval
- **n = 11** : Number of data points
- **p = 0.1** : Point where derivative is to be calculated
- **Function**: f(x) = 1/(1+x²) (hardcoded in code)
- **Step size**: h = (1.0 - 0.0) / 11 ≈ 0.0909

### Output:
```
Difference Table:
      1.0000     -0.0082     -0.0156      0.0022      0.0009     -0.0006      0.0001      0.0001     -0.0001      0.0000      0.0000
      0.9918     -0.0238     -0.0134      0.0031      0.0003     -0.0005      0.0002      0.0000     -0.0001      0.0000
      0.9680     -0.0372     -0.0103      0.0034     -0.0002     -0.0003      0.0002     -0.0000     -0.0000
      0.9308     -0.0476     -0.0069      0.0033     -0.0005     -0.0001      0.0001     -0.0001
      0.8832     -0.0544     -0.0036      0.0028     -0.0007      0.0000      0.0001
      0.8288     -0.0581     -0.0009      0.0021     -0.0007      0.0001
      0.7707     -0.0589      0.0012      0.0014     -0.0006
      0.7118     -0.0577      0.0027      0.0009
      0.6541     -0.0550      0.0035
      0.5990     -0.0515
      0.5475

Forward Interpolation:
1st derivative True value:-0.1961
1st derivative calculated:-0.1963
Error:0.1257%
2nd derivative True value:-1.9018
2nd derivative calculated:-1.8805
Error:1.1176%
```
**Explanation:**
- **Function**: f(x) = 1/(1+x²)
- **True 1st derivative at x=0.1**: f'(0.1) = -2x/(1+x²)² = -0.1961
- **Calculated 1st derivative**: -0.1963 (using forward difference formula)
- **Error**: 0.1257% 
- **True 2nd derivative at x=0.1**: f''(0.1) = -1.9018
- **Calculated 2nd derivative**: -1.8805
- **Error**: 1.1176% 
- **Difference Table**: Shows forward differences up to 10th order
  - First column: function values f(xᵢ)
  - Second column: first differences Δyᵢ
  - Third column: second differences Δ²yᵢ
  - And so on...



---

## Code Constraints

- **Input Format**: Must provide exactly 4 values: b, a, n, p
- **Interval**: **b > a** (upper limit must be greater than lower limit)
- **Data Points**: **n ≥ 5** (recommended for meaningful derivatives, code uses up to 4th order differences)

- **Forward Interpolation**: Works best when **p is close to a** (near beginning of interval)
- **Step Size**: **h = (b-a)/n should be small** for better accuracy
- **Function**: Currently hardcoded as f(x) = 1/(1+x²) in code

- **File Dependency**: Requires input.txt in **../Input/** directory
- **Precision**: Results displayed with 4 decimal places in difference table
- **Memory**: Requires n×n matrix for difference table
- **Order Limitation**: Uses up to 4th order differences; more terms improve accuracy
