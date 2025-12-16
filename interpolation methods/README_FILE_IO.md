# Newton Interpolation with File I/O

## File I/O Feature

The program now supports reading input from `input.txt` and writing output to `output.txt` while simultaneously displaying results in the console.

## Input File Format (input.txt)

```
n y
x[0] f[0]
x[1] f[1]
...
x[n-1] f[n-1]
method_choice
display_choice_1
display_choice_2
...
exit_choice
```

### Example:
```
4 25
20 24
30 30
40 40
50 54
1
1
2
3
```

Where:
- Line 1: `n=4` (number of data points), `y=25` (value to interpolate)
- Lines 2-5: x and f(x) pairs
- Line 6: `1` (Newton Forward) or `2` (Newton Backward)
- Lines 7+: Display choices (`1` = Show Table, `2` = Show Result, `3` = Exit)

## How to Use

1. **Prepare input.txt** with your data and menu choices
2. **Compile**: `g++ -o program main.cpp -std=c++17`
3. **Run**: `./program`
4. **Results**: 
   - Displayed in console (real-time)
   - Saved in `output.txt`

## Files Created

- `file.h` - File I/O handler class
- `NewtonInterpolation.h` - Newton interpolation methods
- `main.cpp` - Main program with file I/O integration

## Features

✓ Reads all input from `input.txt`
✓ Writes output to `output.txt`
✓ Shows results in console simultaneously
✓ No code changes to existing logic
