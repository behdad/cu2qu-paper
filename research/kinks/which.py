import fontTools
from fontTools.cu2qu import cu2qu

print(fontTools.__file__, fontTools.version)
print("cu2qu", cu2qu.__file__, "compiled", cu2qu.COMPILED)
print("has cubic_spline_n", hasattr(cu2qu, "cubic_spline_n"), "QUADRATIC_MAX_D3", getattr(cu2qu, "QUADRATIC_MAX_D3", None))
