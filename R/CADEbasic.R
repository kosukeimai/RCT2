##### basic functions for calculating the CADE

# subtracts element j of list element i of a with element j of list element i of b
# Difflist=function(a,b){
#   J=length(a)
#   c=a
#   for (j in 1:J){
#     c[[j]]=a[[j]]-b[[j]]
#   }
#   return(c)
# }

Meanlist=function(a){
  J=length(a)
  s=0
  for(j in 1:J){
    s=s+mean(a[[j]])
  }
  return(s/J)
}

# multiplies element j of list element i of a with element j of list element i of b
# Productlist=function(a,b){
#   J=length(a)
#   c=a
#   for (j in 1:J){
#     c[[j]]=a[[j]]*b[[j]]
#   }
#   return(c)
# }


Productlist <- function(a, b){
  return( mapply('*', a, b, SIMPLIFY = FALSE) )
}

Difflist <- function(a, b){
  return( mapply('-', a, b, SIMPLIFY = FALSE) )
}