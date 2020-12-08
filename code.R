#2. a)
F<- function(v) v[1]^2 + v[2]^4 - 2*v[3]^2 -2 

v0 <- c(3,1,2)

F(v0)
DF <- grad(F, v0); DF
#Write DF = [A | B] where A = 6 and B = [4, -8]
#The derivative of g is -A^{-1)B
Dg <- -c(DF[2], DF[3])/DF[1]; Dg 

#b)#Suppose we change y to 1.1, z to 1.8. What x satisfies the constraint?
h <- c(0.1, -0.2)    #the increment to y and z
v <- v0 + c(Dg%*%h, h); v
#x is approximately 2.666667

#c) original function x^2+y^4-2*z^2-2 = 0 
#solving for x in terms of y and z 
# g is x = sqrt(-y^4+2*z^2+2)

g <- function(v) sqrt(-v[1]^4 + 2*v[2]^2 + 2)

v0 <- c(1,2)

g(v0)
Dg <- grad(g, v0); Dg

#approximate value of x when y = 1.1 and z = 1.8 
v1 <- c(1.1,1.8)

g(v1)


#3. a)
F <- function(v) c(v[1]*v[2]*v[3]-64, v[1]*v[2]+v[1]*v[3]+v[2]*v[3]-56)
v0 <- c(8, 2, 4)
DF <- jacobian(F, v0); DF
A <- DF[,1:2]; B<- DF[,3]; A; B
AInv <- solve(A); AInv
-AInv%*%B # Dg at the point g=4

#g(4) = c(8,2)
#now calculate approximation for g(3.9) 

c(8,2) + (-AInv%*%B*(3.9-4))


#b. 

F <- function(v) c(v[1]*v[2]*v[3]-64, v[1]*v[2]+v[1]*v[3]+v[2]*v[3]-136)
v0 <- c(2, 2, 16)
DF <- jacobian(F, v0); DF
#specify y, z in terms of x 
A <- DF[,2:3]; B<- DF[,1]; A; B
AInv <- solve(A); AInv
-AInv%*%B # Dg at the point g=2

#specify x, y in terms of z 
A <- DF[,1:2]; B<- DF[,3]; A; B
# AInv <- solve(A); AInv
#A is not invertible! fail to specify x and y in terms of z!



#c)
library("pracma")

F <- function(v) c(v[1]*v[2]*v[3]-64, v[1]*v[2]+v[1]*v[3]+v[2]*v[3]-96)
v0 <- c(4,4,4)
DF <- jacobian(F, v0); DF
rref(DF)
# DF is certainly not onto since it has 3 linearly dependent columns
#therefore there's no implicit function



#4. a)
library(numDeriv)

F <- function(v) c(v[1]^2+v[2]^2 + v[3]^2 - v[4]^2 + 2, 3*v[1]+2*v[2]+v[3]-2*v[4]-2)

#The test for a smooth manifold is that the Jacobian matrix is never [0 0]
#Just plot (or inspect) the first element of each partial derivatives
# since the second row of Jacobian is c(3,2,1,-2)
#never = c(0,0,0,0)

D1F_1 <- function(x) c(2*x)
curve(D1F_1(x), from = -1, to = 1)
D2F_1 <- function(y) c(2*y)
curve(D2F_1(y), from = -1, to = 1, xname = "y", add = TRUE)   #is never zero
D3F_1 <- function(z) c(2*z)
curve(D3F_1(z), from = -1, to = 1, xname = "z", add = TRUE) 
D4F_1 <- function(t) c(-2*t)
curve(D4F_1(t), from = -1, to = 1, xname = "t", add = TRUE) 
abline (h = 0, col = "green")

#all = zero with v = c(0,0,0,0)
# but x^2+y^2+z^2-t^2+2 = 0 
# so the jacobian is onto as the first row is not a scalar multiple of the second row 
#so jacobian is onto and therefore manifold is smooth



# 4. b) 

F <- function(v) c(v[1]^2+v[2]^2 + v[3]^2 - v[4]^2 + 2, 3*v[1]+2*v[2]+v[3]-2*v[4]-2)
v0 <- c(1, 2, 3, 4); F(v0)     #we have a solution
#Near v0, F defines passive variable x implicitly as a function x, y = g(z, t)
#Use the derivative recipe.
DF <- jacobian(F, v0); DF
#Write DF = [A | B]
#The derivative of g is -A^{-1)B
A <- DF[,1:2]; B<- DF[,3:4]; A; B
AInv <- solve(A)
-AInv%*%B # Dg at the point z =3, t=4


#5. 
f<-function(t) 3*t + 2*exp(t) + exp(2*t) -10
df<-function(t) 3+ 2*exp(t) + 2*exp(2*t)

t0 <- 1

t1 <- t0 - f(t0)/df(t0); t1
t2 <- t1 - f(t1)/df(t1); t2
t3 <- t2 - f(t2)/df(t2); t3
t4 <- t3 - f(t3)/df(t3); t4
t5 <- t4 - f(t4)/df(t4); t5

g<-function(t) c(t+exp(t), t+exp(2*t))
g(t5)


library(numDeriv)
# 7. 

f <- function(x,y) x^3 - 12*x*y + 8*y^3   #needed for contour lines
fVec <- function(v) v[1]^3-12*v[1]*v[2] + 8*v[2]^3 #needed for numDeriv
#We can search for maxima and minima by plotting contour lines.
y <- x <- seq(-4, 4, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z) #colors may help identify maximum/minimum.
contour(x,y,z)

#The closed contour for value 0 is suggestive. It is centered near (2,1).
x <- seq(0, 3, .2)
y <- seq(0, 2, 0.2)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)

#It looks as though we have found a minimum near (2, 1).
x <- seq(1.7, 2.5, 0.02)
y <- seq(.8, 1.3, 0.02)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)
contour(x,y,z, levels = c(-7.50,-7.52, -7.54, -7.56, -7.58, -7.59, -7.51))

#Now we have a good approximation.
#Use the partial derivatives to get two equations in two unknowns
Df <- function(v) c(3*v[1]^2-12*v[2], -12*v[1]+24*v[2]^2)
v0 <- c(2, 1)
Df(v0)       #pretty close to zero
grad(fVec, v0)   #in fact, we don't need to know our Math 1 stuff!


H <- hessian(fVec, v0); H; det(H) #positive det confirms extremum
sum(H*diag(c(1,1)))   #positive trace confirms a minimum at c(2,1)!



#zoom into the area where 0 gets very close to itself. It is centered near (0,0).
x <- seq(-1, 1, .2)
y <- seq(-1, 1, 0.2)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)

#It looks as though we have found sth near (0,0).
x <- seq(-.5, .5, 0.02)
y <- seq(-.5, .5, 0.02)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)
contour(x,y,z, levels = c(0.0,-0.02, -0.04, -0.06, -0.08, -0.09, 0.01, 0.02, 0.04))

#Now we have a good approximation.
#Use the partial derivatives to get two equations in two unknowns
Df <- function(v) c(3*v[1]^2-12*v[2], -12*v[1]+24*v[2]^2)
v0 <- c(0, 0)
Df(v0)       #pretty close to zero
grad(fVec, v0)   #in fact, we don't need to know our Math 1 stuff! confirms critical point


H <- hessian(fVec, v0); H; det(H) #negative det, saddle point

#8.
f <- function(x,y) x^3 - y^3 + 2*x*y - 5*x + 6*y   #needed for contour lines
fVec <- function(v) v[1]^3- v[2]^3 + 2*v[1]*v[2] - 5*v[1] + 6*v[2] #needed for numDeriv
#We can search for maxima and minima by plotting contour lines.
y <- x <- seq(-3, 3, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z) #colors may help identify maximum/minimum.
contour(x,y,z)

#critical point #1 

#The closed contour for value 5 is suggestive. It is centered near (-1,1).
x <- seq(-2, 0, .2)
y <- seq(0, 2, 0.2)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)

#It looks as though we have found a minimum near (-1, 1).
x <- seq(-1.2, -.8, 0.02)
y <- seq(1, 1.3, 0.02)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)
contour(x,y,z, levels = c(7,7.02, 7.04, 7.06, 7.08, 7.09, -7.1))

#Now we have a good approximation.
#Use the partial derivatives to get two equations in two unknowns
Df <- function(v) c(3*v[1]^2+2*v[2]-5, -3*v[2]^2+2*v[1]+6)
v0 <- c(-.94, 1.17)
Df(v0)       #pretty close to zero
grad(fVec, v0)   #in fact, we don't need to know our Math 1 stuff!

#Use Newton's method to improve the approximation
A <- jacobian(Df, v0); A
hessian(fVec,v0)  #this matrix of second derivatives can be found numerically
v1 <- v0 - solve(A)%*%Df(v0); v1; Df(v1)
f(v1[1], v1[2])   #the local minimum value of the function

#Look at the Hessian matrix of second partial derivatives at the critical point.
grad(fVec,v1)       #confirms our critical point
H <- hessian(fVec, v1); H; det(H) #positive det confirms extremum
sum(H*diag(c(1,1)))   #negative trace confirms a maximum

#We have found a local minimum. Go back and look for another extremum.
#Replot the original contour lines
y <- x <- seq(-3, 3, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)


#The closed contour for value -10 is suggestive. It is centered near (2, -2).
x <- seq(1, 2, 0.2)
y <- seq(-3, -1, 0.2)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)

#It looks as though we have found a minimum near (1.7, -1.7).
x <- seq(1.6, 1.7, 0.02)
y <- seq(-1.8, -1.7, 0.02)
z <- outer(x,y,f)    #matrix of function values
contour(x,y,z)
contour(x,y,z, levels = c(-14.6,-14.62, -14.64, -14.66, -14.68, -14.69, -14.61))

Df
v0 <- c(1.68, -1.76)
Df(v0)       #pretty close to zero
grad(fVec, v0)   #in fact, we don't need to know our Math 1 stuff!

#Use Newton's method to improve the approximation
A <- jacobian(Df, v0); A
hessian(fVec,v0)  #this matrix of second derivatives can be found numerically
v1 <- v0 - solve(A)%*%Df(v0); v1; Df(v1)
f(v1[1], v1[2])   #the local minimum value of the function

grad(fVec,v1)       #confirms our critical point
H <- hessian(fVec, v1); H; det(H) #positive det confirms extremum
sum(H*diag(c(1,1)))   #positive trace confirms a minimum

#We have found a local max. Go back and look for another extremum.
#Replot the original contour lines
y <- x <- seq(-3, 3, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)

#The 5 contour gets close to itself near (.8,1.5)
abline(h=1.5, col = "green"); abline(v=.8, col = "green")

#Zoom in for a closer look
x <- seq(0, 1.5, 0.2)
y <- seq(1.3, 2, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)

#Something is happening for values between 4.4 and 4.6
contour(x,y,z, levels = seq(4.4, 4.6, 0.01))
#we have found a saddle point near (.8, 1.65)
#Zoom in for a closer look.
x <- seq(.7, .95, 0.01)
y <- seq(1.5, 1.7, 0.01)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)
#Something is happening for values between 4.56 and 4.58.
contour(x,y,z, levels = seq(4.56, 4.58, 0.0005))

#It looks as though we have found a saddle point near (.77, 1.58)
v0 <- c(.77, 1.58)
Df(v0)       #pretty close to zero.
A <- jacobian(Df, v0); A
v1 <- v0 - solve(A)%*%Df(v0); v1; Df(v1) #one iteration of Newton's method works wonders
f(v1[1], v1[2])      #this is the value of the function at the saddle point.
#Notice that there are larger and smaller function values nearby.
#This is not an extremum!

#Have a look at the Hessian matrix 
grad(fVec,v1)       #confirms our critical point
H <- hessian(fVec, v1); H; det(H) #negative determinant confirms saddle point.


#We have found a saddle point. Go back and look for another critical point.
#Replot the original contour lines
y <- x <- seq(-3, 3, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)

#The 0 contour gets close to itself near (-1.5, -1)
abline(h=-1, col = "green"); abline(v=-1.5, col = "green")

#Zoom in for a closer look
x <- seq(-2, -1, 0.2)
y <- seq(-2, 0, 0.2)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)

#Something is happening for values between 2, 2.5
contour(x,y,z, levels = seq(2, 2.5, 0.01))
#we have found a saddle point near (-1.5, -1)
#Zoom in for a closer look.
x <- seq(-1.8, -1.2, 0.01)
y <- seq(-1.5, -.5, 0.01)
z <- outer(x,y,f)    #matrix of function values
filled.contour(x,y,z)
contour(x,y,z)
#Something is happening for values between 2 and 2.2
contour(x,y,z, levels = seq(2, 2.2, 0.0005))

#It looks as though we have found a saddle point near (-1.53, -1)
v0 <- c(-1.53, -1)
Df(v0)       #pretty close to zero.
A <- jacobian(Df, v0); A
v1 <- v0 - solve(A)%*%Df(v0); v1; Df(v1) #one iteration of Newton's method works wonders
f(v1[1], v1[2])      #this is the value of the function at the saddle point.
#Notice that there are larger and smaller function values nearby.
#This is not an extremum!

#Have a look at the Hessian matrix
grad(fVec,v1)       #confirms our critical point
H <- hessian(fVec, v1); H; det(H) #negative determinant confirms saddle point.


