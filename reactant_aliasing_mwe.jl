using Reactant

# Rebinding a local alias should not change the function argument. Under Reactant,
# however, the loop-carried value can alias the input unless it is copied first.
function advance_without_copy(x)
    y = x
    n = 0

    Reactant.@trace track_numbers = true while n < 3
        y += one(y)
        n += 1
    end

    return y
end

function advance_with_copy(x)
    y = copy(x)
    n = 0

    Reactant.@trace track_numbers = true while n < 3
        y += one(y)
        n += 1
    end

    return y
end

# Reuse x after calling the loop function so an unexpected mutation is visible.
without_copy(x) = advance_without_copy(x) + 2x
with_copy(x) = advance_with_copy(x) + 2x

x = [0.25]
expected = 3x .+ 3

bad_input = Reactant.to_rarray(x)
good_input = Reactant.to_rarray(x)

bad = Reactant.@jit without_copy.(bad_input)
good = Reactant.@jit with_copy.(good_input)

println("expected:     ", expected)
println("without copy: ", Array(bad))
println("with copy:    ", Array(good))

@assert Array(good) == expected
@assert Array(bad) != expected
