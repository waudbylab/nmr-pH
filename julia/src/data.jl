# Data loading and selection functions

"""
    load_indicator_data(filename) -> DataFrame

Load chemical shift data from CSV file.

Expected columns:
    - expt_number: experiment identifier (ignored in analysis)
    - ionic_strength: ionic strength in M
    - temperature: temperature in K
    - nominal_pH: pH measured by electrode at 298 K
    - P_xxx, F_xxx, H_xxx: chemical shifts for each indicator
      (prefix indicates nucleus: P=³¹P, F=¹⁹F, H=¹H)

Returns a DataFrame with the data reorganised into long format.
"""
function load_indicator_data(filename)
    raw = CSV.read(filename, DataFrame)

    # Identify indicator columns (those starting with P_, F_, H_)
    indicator_cols = filter(names(raw)) do col
        startswith(string(col), "P_") ||
            startswith(string(col), "F_") ||
            startswith(string(col), "H_")
    end

    # Build long-format dataframe
    rows = []
    for row in eachrow(raw)
        for col in indicator_cols
            col_str = string(col)
            # Extract nucleus and indicator name
            nucleus = col_str[1:1]  # P, F, or H
            indicator = col_str[3:end]  # everything after "X_"

            # Get the chemical shift value
            δ = row[col]

            # Skip missing values
            if δ isa AbstractString
                δ == "missing" && continue
                δ = parse(Float64, string(δ))
            end

            push!(rows, (
                expt_number=row.expt_number,
                ionic_strength=Float64(row.ionic_strength),
                temperature=Float64(row.temperature),
                nominal_pH=Float64(row.nominal_pH),
                nucleus=nucleus,
                indicator=indicator,
                delta_obs=Float64(δ)
            ))
        end
    end

    return DataFrame(rows)
end

"""
    select_condition(data, temperature, ionic_strength; tol_T=1.0, tol_I=0.01) -> DataFrame

Select data for a single (T, I) condition.

Arguments:
    data: DataFrame from load_indicator_data
    temperature: target temperature in K (will select closest)
    ionic_strength: target ionic strength in M (will select closest)
    tol_T: tolerance for temperature matching (K), warns if closest is further
    tol_I: tolerance for ionic strength matching (M), warns if closest is further
"""
function select_condition(data, temperature, ionic_strength; tol_T=1.0, tol_I=0.01)
    # Find closest ionic strength
    I_values = unique(data.ionic_strength)
    closest_I = I_values[argmin(abs.(I_values .- ionic_strength))]

    if abs(closest_I - ionic_strength) > tol_I
        @warn "Requested I=$ionic_strength M, using closest available I=$closest_I M"
    end

    subset = filter(row -> row.ionic_strength == closest_I, data)

    # Find closest temperature
    temps = unique(subset.temperature)
    closest_T = temps[argmin(abs.(temps .- temperature))]

    if abs(closest_T - temperature) > tol_T
        @warn "Requested T=$temperature K, using closest available T=$closest_T K"
    end

    subset = filter(row -> row.temperature == closest_T, subset)

    println("Selected $(nrow(subset)) measurements at T=$(closest_T) K, I=$(closest_I) M")
    println("Indicators: ", join(sort(unique(subset.indicator)), ", "))

    return subset
end

"""
    list_conditions(data)

Print available (temperature, ionic_strength) combinations in the data.
"""
function list_conditions(data)
    conditions = unique(select(data, [:temperature, :ionic_strength]))
    sort!(conditions, [:temperature, :ionic_strength])

    println("\nAvailable conditions:")
    println("─"^30)
    for row in eachrow(conditions)
        n = count(r -> r.temperature == row.temperature &&
                r.ionic_strength == row.ionic_strength, eachrow(data))
        println("  T = $(row.temperature) K, I = $(row.ionic_strength)  ($n measurements)")
    end
    println()
end

"""
    get_indicator_data(data, indicator) -> (pH_values, delta_values, nucleus)

Extract pH and chemical shift values for a single indicator.
"""
function get_indicator_data(data, indicator::String)
    subset = filter(row -> row.indicator == indicator, data)
    if nrow(subset) == 0
        error("Indicator '$indicator' not found in data")
    end
    nucleus = subset.nucleus[1]
    pH_vals = Float64.(subset.nominal_pH)
    delta_vals = Float64.(subset.delta_obs)
    return pH_vals, delta_vals, nucleus
end
