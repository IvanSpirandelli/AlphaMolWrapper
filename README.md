This is a Julia wrapper for [AlphaMol](https://github.com/pkoehl/AlphaMol) from [Patrice Koehl](https://www.cs.ucdavis.edu/~koehl/)

Its purpose is to calculate intrinsic volumes of a union of balls in three dimensions.

Example of calling one of the wrapped functions from Julias side:

```
module AlphaMolWrap
    using CxxWrap
    using AlphaMolWrapper_jll

    @wrapmodule(() -> libalphamolwrapper)

    function __init__()
        @initcxx
    end
end

function get_measures(coordinates::Vector, radii::Vector; delaunay_eps::Float64 = 1.0)
    measures_out = [0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.calculate_measures(
        measures_out,
        coordinates,
        radii,
        delaunay_eps, 
        0, # Calculate derivatives flag. 0 for no 1 for yes
        0  # Print debug output flag. 0 for no 1 for yes 
    )
    measures_out
end
```