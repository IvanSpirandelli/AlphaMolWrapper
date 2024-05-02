This is a Julia wrapper for [AlphaMol](https://github.com/pkoehl/AlphaMol) from [Patrice Koehl](https://www.cs.ucdavis.edu/~koehl/)

Its purpose is to calculate intrinsic volumes of a union of balls in three dimensions.

We added a calculation of overlap values that allows for the computation of overlaps between different molecules of the same type or individual atoms. 
It utilizes the existing Alpha Complex to check for violated boundaries and sums up overlap penalties.

Examples of calling the wrapped functionality from Julias side:

```
module AlphaMolWrap
    using CxxWrap
    using AlphaMolWrapper_jll

    @wrapmodule(() -> libalphamolwrapper)

    function __init__()
        @initcxx
    end
end

function get_geometric_measures_and_overlap_value(
    atom_coordinates::Vector, 
    molecule_size::Int, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64, 
    delaunay_eps::Float64 = 1.0)
    outs = [0.0, 0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.get_geometric_measures_and_overlap_value(
        outs,
        atom_coordinates,
        molecule_size,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    outs
end

function get_geometric_measures_and_overlap_value_with_derivatives(
    atom_coordinates::Vector, 
    molecule_size::Int, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    overlap_existence_penalty::Float64,
    overlap_penalty_slope::Float64, 
    delaunay_eps::Float64 = 1.0)
    
    measure_outs = [0.0, 0.0, 0.0, 0.0, 0.0]

    n = size(atom_coordinates)[1]

    dvol_outs = [0.0 for _ in 1:n]
    dsurf_outs = [0.0 for _ in 1:n]
    dmean_outs = [0.0 for _ in 1:n]
    dgauss_outs = [0.0 for _ in 1:n]
    dlol_outs = [0.0 for _ in 1:n]

    AlphaMolWrap.get_geometric_measures_and_overlap_value_with_derivatives(
        measure_outs,
        dvol_outs,
        dsurf_outs, 
        dmean_outs,
        dgauss_outs,
        dlol_outs,
        atom_coordinates,
        molecule_size,
        atom_radii,
        probe_radius,
        overlap_existence_penalty,
        overlap_penalty_slope,
        delaunay_eps
    )
    measure_outs, dvol_outs, dsurf_outs, dmean_outs, dgauss_outs, dlol_outs
end

function get_geometric_measures(
    atom_coordinates::Vector, 
    atom_radii::Vector, 
    probe_radius::Float64,
    delaunay_eps::Float64 = 1.0)
    outs = [0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.get_geometric_measures(
        outs,
        atom_coordinates,
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    outs
end

function get_geometric_measures_with_derivatives(
    atom_coordinates::Vector, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    delaunay_eps::Float64 = 1.0)
    
    measure_outs = [0.0, 0.0, 0.0, 0.0, 0.0]

    n = size(atom_coordinates)[1]

    dvol_outs = [0.0 for _ in 1:n]
    dsurf_outs = [0.0 for _ in 1:n]
    dmean_outs = [0.0 for _ in 1:n]
    dgauss_outs = [0.0 for _ in 1:n]

    AlphaMolWrap.get_geometric_measures_with_derivatives(
        measure_outs,
        dvol_outs,
        dsurf_outs, 
        dmean_outs,
        dgauss_outs,
        atom_coordinates,
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    measure_outs, dvol_outs, dsurf_outs, dmean_outs, dgauss_outs
end

```

The parameter `delaunay_eps` defines the precision with which delaunay triangulations are calculated when computing the geometric measures. Lowering this value will speed up calculations but might reduce accuracy. The perfect value depends on the scale of the calculated unions. 