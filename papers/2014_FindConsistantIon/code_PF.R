
#data <- read.csv("F:data.csv", row.names=1)

data <- read.table("F:data.txt", header = TRUE)

#> dim(data)
#[1] 50557 12470


##Calculate p-value 

z.pvalue <- function( x ){
              se <- sd( x )/sqrt( length( x ) )
              xbar <- mean( x )
              z <- ( xbar - 0 )/se
              pnorm( z, mean = 0, sd = 1, lower.tail = FALSE)*2
}


data1 <- data[1:50557,]

p_table <- data.frame( apply( data1, 1, z.pvalue ) )
colnames(p_table) <-"p_value"
p_table$ion <- c(data1$Ion)
p_table$adj_p <- p.adjust(p_table$p_value)
Final_table <- data.frame(p_table[ order( p_table$p_value),])

AA <- head(Final_table[,c("ion","p_value","adj_p")], 5000)

# write.csv(AA, file = "F:ION_P_value.csv")

