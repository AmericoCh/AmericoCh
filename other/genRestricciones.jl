# Algoritmo 1: Generación de restricciones para RO. 
# Resolvemos problemas determinísticos
# con un número creciente de restricciones 
# eficiente?
using JuMP
using Gurobi
function robust_optimization(f, g, ξ0, U)
    # Inicialización
    k = 1
    ξ_list = [ξ0]
    x_opt = nothing
    v_approx = nothing

    # Iteración principal
    while true
        # Resolver el subproblema restringido (optimización en x)
        model_min = Model(Gurobi.Optimizer)
        @variable(model_min, x[1:2])  
        @objective(model_min, Min, f(x))

        for ξ in ξ_list
            @constraint(model_min, g(x, ξ) <= 0)
        end
        optimize!(model_min)

        if termination_status(model_min) != MOI.OPTIMAL
            error("No se pudo encontrar una solución óptima.")
        end
        x_k = value.(x)
        v_approx = objective_value(model_min)

        # Resolver la maxim de la violación de restricciones (optimización sobre ξ)
        model_max = Model(Gurobi.Optimizer)
        @variable(model_max, ξ[1:length(ξ0)])
        for i in 1:length(ξ0)
            @constraint(model_max, U[i][1] <= ξ[i] <= U[i][2]) 
        end
        @objective(model_max, Max, g(x_k, ξ))
        optimize!(model_max)

        if termination_status(model_max) != MOI.OPTIMAL
            error("No se pudo encontrar una solución óptima.")
        end

        s = objective_value(model_max)
        ξ_k_new = value.(ξ)

        if s <= 0
            return v_approx, x_k
        end

        push!(ξ_list, ξ_k_new)
        k += 1
    end
end

# Ejemplo 
f(x) = x[1]^2 + x[2]^2  # Función objetivo
g(x, ξ) = ξ[1] * x[1] + ξ[2] * x[2] - 1  # Restricción robusta
ξ0 = [0.5, 0.5]
U = [[-1.0, 1.0], [-1.0, 1.0]]  # Conjunto de incertibumbre

v, x = robust_optimization(f, g, ξ0, U)
println("Valor aproximado: $v, Solución: $x")
