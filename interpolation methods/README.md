# Newton's Interpolation Methods

## Theory

Newton's Interpolation is a method of constructing an interpolating polynomial that passes through a given set of data points. It uses the concept of divided differences (for unequal intervals) or finite differences (for equal intervals) to estimate values between known data points.

### Mathematical Basis:

Interpolation is used to find intermediate values of a function when only discrete data points are known. Newton developed two main approaches:

1. **Forward Interpolation**: Used when x is near the beginning of the data table
2. **Backward Interpolation**: Used when x is near the end of the data table

---

## Newton's Forward Interpolation

### When to Use:
- Forward interpolation is used when **x is near the beginning** of the table
- Best for interpolating values in the first half of the data range

### For Equal Intervals:

Let x₀, x₁, x₂, ..., xₙ₋₁, xₙ be a set of equidistant values of variable x.

**Equal spacing**: x₁ - x₀ = x₂ - x₁ = x₃ - x₂ = ... = xₙ - xₙ₋₁ = h

**Let**: u = (x - x₀) / h

### Newton's Forward Difference Formula:

**For Equal Intervals**:
```
y = y₀ + uΔy₀ + [u(u-1)/2!]Δ²y₀ + [u(u-1)(u-2)/3!]Δ³y₀ + ... + [u(u-1)(u-2)...(u-n+1)/n!]Δⁿy₀
```

Where:
- Δy₀ = y₁ - y₀ (first forward difference)
- Δ²y₀ = Δy₁ - Δy₀ (second forward difference)
- Δⁿy₀ = nth forward difference

**For Unequal Intervals** (General Formula):
```
f(xₙ) = f(x₀) + (x-x₀)f[x₁,x₀] + (x-x₀)(x-x₁)f[x₂,x₁,x₀] + ... 
        + (x-x₀)(x-x₁)...(x-xₙ₋₁)f[xₙ,xₙ₋₁,...,x₁,x₀]
```

Where divided differences are:
```
f[xᵢ,xⱼ] = [f(xᵢ) - f(xⱼ)] / (xᵢ - xⱼ)

f[xᵢ,xⱼ,xₖ] = [f[xᵢ,xⱼ] - f[xⱼ,xₖ]] / (xᵢ - xₖ)
```

### Error Calculation:
If an additional data point f(xₙ₊₁) is available, the truncation error can be approximated as:
```
Rₙ ≈ f[xₙ₊₁,xₙ,xₙ₋₁,...,x₁,x₀](x - x₀)(x - x₁)...(x - xₙ)
```

---

## Newton's Backward Interpolation

### When to Use:
- Backward interpolation is used when **x is near the end** of the table
- Best for interpolating values in the last half of the data range

### For Equal Intervals:

Let x₀, x₁, x₂, ..., xₙ₋₁, xₙ be a set of equidistant values of variable x.

**Equal spacing**: x₁ - x₀ = x₂ - x₁ = x₃ - x₂ = ... = xₙ - xₙ₋₁ = h

**Let**: v = (x - xₙ) / h

### Newton's Backward Difference Formula:
```
y = yₙ + v∇yₙ + [v(v+1)/2!]∇²yₙ + [v(v+1)(v+2)/3!]∇³yₙ + ... 
    + [v(v+1)(v+2)...(v+n-1)/n!]∇ⁿyₙ
```

Where:
- ∇yₙ = yₙ - yₙ₋₁ (first backward difference)
- ∇²yₙ = ∇yₙ - ∇yₙ₋₁ (second backward difference)
- ∇ⁿyₙ = nth backward difference

---

## Comparison: Forward vs Backward

| Aspect | Forward Interpolation | Backward Interpolation |
|--------|----------------------|------------------------|
| **Starting Point** | x₀ (beginning) | xₙ (end) |
| **Best Used When** | x near start of table | x near end of table |
| **Difference Operator** | Δ (forward difference) | ∇ (backward difference) |
| **Parameter** | u = (x-x₀)/h | v = (x-xₙ)/h |
| **Table Direction** | Top to bottom | Bottom to top |

### Algorithm:

**Step 1: Input Data**
- Read n data points (xᵢ, yᵢ)
- Read the value x at which to interpolate

**Step 2: Determine Method**
- If x is closer to x₀: Use Forward Interpolation
- If x is closer to xₙ: Use Backward Interpolation

**Step 3: Construct Difference Table**
- **Forward**: Calculate Δy, Δ²y, Δ³y, ..., starting from y₀
- **Backward**: Calculate ∇y, ∇²y, ∇³y, ..., starting from yₙ

**Step 4: Calculate Parameter**
- **Forward**: u = (x - x₀) / h
- **Backward**: v = (x - xₙ) / h

**Step 5: Apply Formula**
- Compute y using the appropriate Newton formula
- Sum all terms until desired accuracy

**Step 6: Output Result**
- Display interpolated value f(x)


---

## Input/Output Example

### Input Format:
```
n y
x₀ f(x₀)
x₁ f(x₁)
x₂ f(x₂)
...
xₙ f(xₙ)
```

### Input:
```
4 3.8
1 1
2 8
3 27
4 64
```
**Explanation:**
- n = 4 (4 data points)
- y = 3.8 (value to interpolate at)
- Data points: (1,1), (2,8), (3,27), (4,64)
- Function appears to be f(x) = x³
- y = 3.8 is near the end, so backward interpolation is better

### Output:
```
=== Backward Difference Table ===
1.00 
8.00 7.00 
27.00 19.00 12.00 
64.00 37.00 18.00 6.00 

=== Newton Backward Result ===
f(3.8) = 54.872000
```
**Explanation:**
- **Backward Difference Table**:
  - Column 0: Original f(x) values: 1, 8, 27, 64
  - Column 1: 1st differences: 7, 19, 37
  - Column 2: 2nd differences: 12, 18
  - Column 3: 3rd differences: 6
- **Result**: f(3.8) ≈ 54.872
- **Verification**: 3.8³ = 54.872 
- **Why Backward**: x=3.8 is close to x₃=4 (end of table)

---


## Code Constraints:

- **Input Format**: First line: n y, then n lines of xᵢ f(xᵢ)
- **Equal Spacing Required**: x₁-x₀ = x₂-x₁ = ... = h
- **Data Points**: n ≥ 2 (at least 2 points needed)
- **Sorted Data**: x values should be in ascending order
- **Range**: y should ideally be within [x₀, xₙ] for best accuracy
- **File Dependency**: Uses file I/O for input/output
- **Memory**: O(n²) for difference tables
