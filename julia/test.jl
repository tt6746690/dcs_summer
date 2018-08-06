
function simulate_π(n)
    in_quad_circle(x, y) = x^2 + y^2 <= 1
    τ = 0
    for i = 1:n
        x, y = rand(), rand()
        if in_quad_circle(x, y)
            τ += 1
        end
    end
    π = 4*τ / n
    return π
end

π = simulate_π(10000)
println(π)  # 3.142
