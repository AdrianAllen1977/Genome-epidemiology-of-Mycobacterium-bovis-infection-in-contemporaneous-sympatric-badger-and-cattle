# Create two tables
geneticDistances <- matrix(rnorm(n=100), nrow=10, ncol=10)
hostRelatedness <- matrix(rnorm(n=100), nrow=10, ncol=10)

# Initialise a table to store the paired values of these tables as columns
table <- data.frame("Genetic"=NA, "HostRelatedness"=NA)
row <- 0

# Fill the vector by selecting each value in the lower triangle
for(i in 1:nrow(geneticDistances)){
        for(j in 1:ncol(geneticDistances)){

                # Skip the upper triangle and self comparisons
                if(i >= j){
                        next
                }

                # Print where you currently are in table
                cat(paste("Current row = ", i, "\tCurrent column = ", j, "\n", sep=""$

                # Store current value
                row <- row + 1
                table[row, "Genetic"] <- geneticDistances[i, j]
                table[row, "HostRelatedness"] <- hostRelatedness[i, j]
        }
}
#### Function  -get upper traingale of matrix.
getUpperTriangle <- function(matrix){

  

  # Initialise a vector to store the values in the upper triangle

  vector <- c()

  

  # Use nested loops to visit each entry in upper trianle

  for(i in 1:nrow(matrix)){

    for(j in 1:ncol(matrix)){

      

      # Ignore upper triangle and self comparisons

      if(i >= j){

        next

      }

      

      # Note progress

      #Sys.sleep(1) # Make computer sleep for 1 second

      #cat(paste("Current row =", i, "\tCurrent column =", j, "\n"))

      

      # Store value

      vector[length(vector) + 1] <- matrix[i, j]

    }

  }

  

  return(vector)

}
