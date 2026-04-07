using Catalyst, Plots, OrdinaryDiffEq
essaie(A,x) = A*x/A
@species A(t) B(t) C(t) D(t) Y(t) X(t) 
rn1 = @network_component rn1 begin 
        @species A(t) 


      k1, A --> B
     k2, B --> C

end


rn2 = @network_component rn2 begin
    @species A(t) 

      k3, A --> D  
      k4, C --> 0

end 



rn3 = @network_component rn3 begin  
    @species A(t) 
    essaie(A,k4), Y--> X
    k5, X--> 0

end


CI = [:A => 0 , :B=> 1000,:C=> 0,:D=> 0.0,:X=> 0,:Y=> 200]
print(CI)

p = [:k1 => 0.1 , :k2=>0.1, :k3=>0.01,:k4=>0.2,:k5=>0.002]

rn_full = extend(rn1,rn2)
rn_full_full = extend(rn_full,rn3)
t_span =(0,360)

prob = ODEProblem(complete(rn_full_full),CI,t_span,p)
sol = solve(prob)

plot(sol.t,sol[:D],label="D")