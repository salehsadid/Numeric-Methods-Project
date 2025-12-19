# Simpson's 3/8ᵗʰ Rule

## Theory

Simpson's 3/8 Rule is part of the Newton-Cotes family of numerical integration formulas. It is based on replacing a function with a polynomial that passes through evenly spaced points. Polynomials are used because they are simple to integrate.

**Goal:** Approximate definite integrals when the exact integral is difficult or impossible to compute analytically.

### Classification:
Simpson's 3/8 Rule is a **closed-form Newton-Cotes formula** (endpoints included) that uses **4 points** to fit a **cubic polynomial** (third-degree) through equally spaced data points.

![Simpson's 3/8 Rule Illustration](../../RESOURCES/threeeight.png)

### Key Properties:
- Uses **4 equally spaced points** (requires number of intervals to be a multiple of 3)
- Approximates the curve with a third-degree (cubic) polynomial
- Useful when the total number of intervals is a multiple of 3
- Slightly more accurate than Simpson's 1/3 Rule in some cases

### Formula:
```
∫ₐᵇ y dx = (3h/8) × [(y₀ + yₙ) + 3(y₁ + y₂ + y₄ + ⋯ + yₙ₋₁) + 2(y₃ + y₆ + ⋯ + yₙ₋₃)]
```

Where:
- **h** = interval width = (b - a) / n
- **b** = upper limit of integration
- **a** = lower limit of integration
- **n** = number of equal parts (must be a multiple of 3)
- **y₀, y₁, y₂, ..., yₙ** = function values at equally spaced points


---

## Algorithm

**Step 1:** Read the number of intervals **n** (must be a multiple of 3), lower limit **a**, and upper limit **b**

**Step 2:** Calculate the interval width: **h = (b - a) / n**

**Step 3:** Calculate function values at all points:
- y₀ = f(a)
- y₁ = f(a + h)
- y₂ = f(a + 2h)
- ...
- yₙ = f(b)

**Step 4:** Separate the function values:
- sum₁ = y₀ + yₙ (first and last)
- sum₂ = 3 × (y₁ + y₂ + y₄ + y₅ + y₇ + y₈ + ... + yₙ₋₁) (all non-multiple-of-3 indices)
- sum₃ = 2 × (y₃ + y₆ + y₉ + ... + yₙ₋₃) (multiple-of-3 indices except endpoints)

**Step 5:** Apply Simpson's 3/8 formula:
```
Result = (3h/8) × (sum₁ + sum₂ + sum₃)
```

**Step 6:** Output the result

---

## Input/Output Example

### Input Format:
```
n a b
```

### Input:
```
6 1 4
```
**Explanation:**
- **n = 6** : Number of intervals (must be a multiple of 3)
- **a = 1** : Lower limit of integration
- **b = 4** : Upper limit of integration
- Function integrated: f(x) = √x

### Output:
```
f(x) = sqrt(x)
a = 1, b = 4, n = 6, h = 0.5
Result = 4.666461
```
**Explanation:**
- Interval width h = (4-1)/6 = 0.5
- Points: x₀=1.0, x₁=1.5, x₂=2.0, x₃=2.5, x₄=3.0, x₅=3.5, x₆=4.0
- Applied Simpson's 3/8 formula: ∫₁⁴ √x dx ≈ 4.666461


---

## Code Constraints

- **Input Format**: Must provide exactly 3 values: n, a, b
- **Interval Count**: **n must be a multiple of 3** (required for Simpson's 3/8 Rule)
- **Limits**: **a < b** (lower limit must be less than upper limit)
- **Interval Width**: **h = (b-a)/n must be > 0**
- **File Dependency**: Requires input.txt in **../Input/** directory
- **Function**: Currently hardcoded as f(x) = √x in code

- **Comparison**: Slightly more accurate than Simpson's 1/3 Rule but requires more function evaluations
