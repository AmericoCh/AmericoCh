#Algoritmo de Benders  para modelos ARSO
#Caso de estudio: Expansión de la Transmisión
#07/01/25
# version 1.1  # Automatizar add cuts?!



using JuMP
using Gurobi



const M = 10000  # Constante grande para restricciones de variables binarias ()
const TOLERANCE = 1e-4  # Tolerancia para la convergencia
const MAX_ITERATIONS = 100  # Número máximo de iteraciones



#PROBLEMA MAESTRO
function solve_master_problem(ps_prev, pc_prev, d_prev)
    m = Model(Gurobi.Optimizer)

    
    # Varaibles Binarias
    @variable(m, x1, Bin)
    @variable(m, x2, Bin)
    
    #Variables positivas
    @variable(m, eta >= 0)
    @variable(m, ps11 >= 0)
    @variable(m, pc11 >= 0)
    @variable(m, pu11 >= 0)
    @variable(m, ps21 >= 0)
    @variable(m, pc21 >= 0)
    @variable(m, pu21 >= 0)
    
    
    #Variables irrestrictas
    @variable(m, f111)
    @variable(m, f211)
    @variable(m, d111)
    @variable(m, d211)
    @variable(m, f121)
    @variable(m, f221)
    @variable(m, d121)
    @variable(m, d221)
    

    #Función objetivo
    @objective(m, Min, 2*10^6 * x1 + 3*10^6 * x2 + eta)
    

    # Restricciones
    @constraint(m, eta >= 0.5 * 8760 * (2 * ps11 + 20 * pc11 + 200 * pu11) +
                      0.5 * 8760 * (2 * ps21 + 20 * pc21 + 200 * pu21))
    
    # Escenario 1
    @constraint(m, ps11 == f111)
    @constraint(m, pc11 == f211)
    @constraint(m, 0.5 * d_prev - pu11 == f111 + f211)
    @constraint(m, f111 - 100 * d111 >= -M * (1 - x1))  
    @constraint(m, f111 - 100 * d111 <= M * (1 - x1))   
    @constraint(m, f211 - 100 * d211 >= -M * (1 - x2)) 
    @constraint(m, f211 - 100 * d211 <= M * (1 - x2))  
    @constraint(m, -50*x1 <= f111)
    @constraint(m,f111 <= 50*x1)
    @constraint(m, -50*x2 <= f211)
    @constraint(m,f211 <= 50*x2)
    @constraint(m, ps11 <= 0.75 * ps_prev)
    @constraint(m, pc11 <= 1.00 * pc_prev)
    @constraint(m, pu11 <= 0.50 * d_prev)
    
    # Escenario 2
    @constraint(m, ps21 == f121)
    @constraint(m, pc21 == f221)
    @constraint(m, 0.75 * d_prev - pu21 == f121 + f221)
    @constraint(m, f121 - 100 * d121 >= -M * (1 - x1))  
    @constraint(m, f121 - 100 * d121 <= M * (1 - x1))   
    @constraint(m, f221 - 100 * d221 >= -M * (1 - x2))  
    @constraint(m, f221 - 100 * d221 <= M * (1 - x2))   
    @constraint(m, -50*x1 <= f121)
    @constraint(m, f121 <= 50*x1)
    @constraint(m, -50*x2 <= f221)
    @constraint(m, f221 <= 50*x2)
    
    @constraint(m, ps21 <= 0.25 * ps_prev)
    @constraint(m, pc21 <= 1.00 * pc_prev)
    @constraint(m, pu21 <= 0.75 * d_prev)

    optimize!(m)

    return value(eta), value(x1), value(x2)
end



#SUBPROBLEMA
function solve_subproblem(x1, x2)
    s = Model(Gurobi.Optimizer)
    set_optimizer_attribute(s, "NonConvex", 2) # resolver prob no convexos con flexibilidad

    #variables
    @variable(s, PS)
    @variable(s, PC)
    @variable(s, D)
    
    #variables irrestrictas: duales de restricciones de igualdad en el Escenario 1
    @variable(s, alpha11)
    @variable(s, alpha12)
    @variable(s, alpha13)
    @variable(s, alpha21)
    @variable(s, alpha22)
    @variable(s, alpha23)
    
    
    #Variables positivas: duales  de restricciones de desigualdad en el Escenario 1
    @variable(s, betamax111 >= 0)
    @variable(s, betamin111 >= 0)
    @variable(s, betamax112 >= 0)
    @variable(s, betamin112 >= 0)
    @variable(s, betamax121 >= 0)
    @variable(s, betamin121 >= 0)
    @variable(s, betamax122 >= 0)
    @variable(s, betamin122 >= 0)
    @variable(s, gammamax11 >= 0)
    @variable(s, gammamin11 >= 0)
    @variable(s, gammamax12 >= 0)
    @variable(s, gammamin12 >= 0)
    @variable(s, gammamax13 >= 0)
    @variable(s, gammamin13 >= 0)
    
    
    #variables positivas: duales de restricciones de desigualdad en el Escenario 2
    @variable(s,betamax211, lower_bound = 0)
    @variable(s,betamin211, lower_bound =0)
    @variable(s,betamax212, lower_bound =0)
    @variable(s,betamin212, lower_bound =0)
    @variable(s,betamax221, lower_bound =0)
    @variable(s,betamin221, lower_bound =0)
    @variable(s,betamax222, lower_bound =0)
    @variable(s,betamin222, lower_bound =0)
    @variable(s,gammamax21, lower_bound =0)
    @variable(s,gammamin21, lower_bound =0)
    @variable(s,gammamax22, lower_bound =0)
    @variable(s,gammamin22, lower_bound =0)
    @variable(s,gammamax23, lower_bound =0)
    @variable(s,gammamin23, lower_bound =0)
      
      
    
     #Funcion objetivo
    @objective(s, Max,-( (1-x1)*10e6*betamax111 + (1-x1)*10e6*betamin111 + x1*50*betamax112 + x1*50*betamin112 
                        + (1-x2)*10e6*betamax121 + (1-x2)*10e6*betamin121 + x2*50*betamax122 + x2*50*betamin122
                        + 0.75*PS*gammamax11 +1.0* PC*gammamax12 + 0.5*D*gammamax13 - 0.5*D*alpha13
                        + (1-x1)*10e6* betamax211 + (1-x1)* 10e6* betamin211 + x1*50*betamax212 + x1*50*betamin212
                        + (1-x2)*10e6*betamax221 + (1-x2)*10e6*betamin221 + x2*50*betamax222 + x2*50*betamin222
                        + 0.25*PS*gammamax21 + 1.0*PC*gammamax22 + 0.75*D*gammamax23 - 0.75*D*alpha23))

    
    # Restricciones de variables inciertas (incertibumbre a largo plazo, problema de segundo nivel) 
    @constraint(s, 140 <= PS <= 200)
    @constraint(s, 80 <= PC <= 100)
    @constraint(s, 64 <= D <= 80) 
    @constraint(s,200 - PS <= 60)   #Modificar el presupuesto de incertidumbre?!
    @constraint(s,100 - PC <= 20)
    @constraint(s, D - 64 <= 16)
    
    
    #Restricciones duales del Escenario 1
    @constraint(s, 0.5*8760*2 + alpha11 + gammamax11 - gammamin11 == 0)
    @constraint(s, 0.5*8760*20 + alpha12 + gammamax12 - gammamin12 == 0)
    @constraint(s, 0.5*8760*200 - alpha13 + gammamax13 - gammamin13 == 0)

    @constraint(s, -alpha11 - alpha13 - betamin111/100 +betamax111/100 - betamin112 + betamax112 == 0)
    @constraint(s, -alpha12 - alpha13 - betamin121/100 +betamax121/100 - betamin122 + betamax122 == 0)

    @constraint(s, betamin111 - betamax111 == 0)
    @constraint(s, betamin121 - betamax121 == 0)



    #Restricciones duales del segundo Escenario
    @constraint(s, 0.5*8760*2 + alpha21 + gammamax21 - gammamin21 == 0)
    @constraint(s, 0.5*8760*20 + alpha22 + gammamax22 - gammamin22 == 0)
    @constraint(s, 0.5*8760*200 - alpha23 + gammamax23 - gammamin23 == 0)

    @constraint(s, -alpha21 - alpha23 - betamin211/100 +betamax211/100 - betamin212 + betamax212 == 0)
    @constraint(s, -alpha22 - alpha23 - betamin221/100 +betamax221/100 - betamin222 + betamax222 == 0)

    @constraint(s, betamin211 - betamax211 == 0)
    @constraint(s, betamin221 - betamax221 == 0)
    
    

    optimize!(s)

    return value(PS), value(PC), value(D)
end




function benders_algorithm()
    k = 1
    eta_prev = 1e9  # Cota superior inicial
    eta_lb = 0  # Cota inferior inicial
    ps, pc, d = 150.0, 90.0, 70.0  # Valores iniciales para variables inciertas 6 2 2???

    while k <= MAX_ITERATIONS && abs(eta_prev - eta_lb) > TOLERANCE
        println("Iteración $k:")
        eta, x1, x2 = solve_master_problem(ps, pc, d)
        ps, pc, d = solve_subproblem(x1, x2)

        # Actualizar cotas
        eta_lb = eta
        eta_prev = min(eta_prev, eta)

        println("Cota inferior: $eta_lb, Cota superior: $eta_prev")
        println("x1: $x1, x2: $x2")
        k += 1
    end
end

# Ejecutar el algoritmo
benders_algorithm()
