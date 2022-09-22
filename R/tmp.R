

# -----
# compute cluster size parameter
# Use PCA to compute eigenvector of geometry this provides ~diamter of cluster




# ----
# Object found using localizations and not merged localizations tend to be bigger!

p <- ggplot(data = cls_counts[cls_counts$n>3, ], aes(x = n, y = PCA1_nm))
p <- p + geom_point()
p

# ----

# Plot one of the clusters in 3D


data <-read.csv(csv_files[1])
cls_counts <- data %>% count(data[cls])
# find biggest cluster
big_cls <- cls_counts[which.max(cls_counts$n), cls]

# PCA
dd <- dist(data[data[cls]==big_cls, c('ltr_x', 'ltr_y', 'ltr_z') ])
max_size <- max(dd)
PCA <- prcomp(data[data[cls]==big_cls, c('ltr_x', 'ltr_y', 'ltr_z') ])
x<-list()
y <- list()
z <- list()
for (id_vec in c(1,2,3)){
  
  x <- rbind(x, cbind(as.vector(PCA$center[1] - PCA$rotation[1,id_vec]*PCA$sdev[id_vec]), as.vector(PCA$center[1] + PCA$rotation[1,id_vec]*PCA$sdev[id_vec])))
  y <- rbind(y, cbind(as.vector(PCA$center[2] - PCA$rotation[2,id_vec]*PCA$sdev[id_vec]), as.vector(PCA$center[2] + PCA$rotation[2, id_vec]*PCA$sdev[id_vec])))
  z <- rbind(z, cbind(as.vector(PCA$center[3] - PCA$rotation[3,id_vec]*PCA$sdev[id_vec]), as.vector(PCA$center[3] + PCA$rotation[3, id_vec]*PCA$sdev[id_vec])))
}

x <- rbind(x, cbind(as.vector(PCA$center[1] - PCA$rotation[1, 1]*max_size/2), 
                    as.vector(PCA$center[1] + PCA$rotation[1, 1]*max_size/2)))
y <- rbind(y, cbind(as.vector(PCA$center[2] - PCA$rotation[2, 1]*max_size/2), 
                    as.vector(PCA$center[2] + PCA$rotation[2, 1]*max_size/2)))
z <- rbind(z, cbind(as.vector(PCA$center[3] - PCA$rotation[3, 1]*max_size/2), 
                    as.vector(PCA$center[3] + PCA$rotation[3, 1]*max_size/2)))

x 
scene <- list(camera = list(eye = list(x = -1.25, y = 1.25, z = 1.0)), 
              xaxis = list(title='x'), 
              yaxis = list(title = 'y'), 
              zaxis=list(title = 'z'))

fig <- plot_ly() %>% add_trace(data=data[data[cls]==big_cls, ], x = ~ltr_x, y = ~ltr_y, z = ~ltr_z, #, size=~se, sizes=c(5, 50), 
                               marker = list(symbol = 'circle', sizemode = 'diameter'), 
                               mode='markers', type='scatter3d',  name=paste('Cluster', big_cls) )

fig <- fig %>% add_trace(x=~x[1,] , y=~y[1,], z=~z[1,],  type="scatter3d", mode="lines", line = list(width=10), 
                         showlegend = TRUE,  name=paste('Longest axis', round(2*PCA$sdev[1]*1e9, 1), '(nm)')) 
fig <- fig %>% add_trace(x=~x[2,] , y=~y[2,], z=~z[2,], type="scatter3d", mode="lines", line = list(width=10), showlegend = FALSE) 
fig <- fig %>% add_trace(x=~x[3,] , y=~y[3,], z=~z[3,], type="scatter3d", mode="lines", line = list(width=10), showlegend = FALSE) 
fig <- fig %>% add_trace(x=~x[4,] , y=~y[4,], z=~z[4,], type="scatter3d", mode="lines", line = list(width=10), showlegend = FALSE) 

fig <- fig %>% layout(scene=scene)
fig
#saveWidget(fig, "C:/Users/apoliti/Desktop/temp.html")
#webshot("C:/Users/apoliti/Desktop/temp.html", "C:/Users/apoliti/Desktop/temp.png")
# ----
fig <- plot_ly(data=data[data[cls]==big_cls, ], x = ~ltr_x, y = ~ltr_y, z = ~ltr_z,  
               marker = list(symbol = 'circle', sizemode = 'diameter'), 
               type='scatter3d',  name=paste('Cluster', big_cls) )
fig <- fig %>% add_markers()
#fig <- fig %>% add_trace(x=~x[1,] , y=~y[1,], z=~z[1,],  type="scatter3d", mode="lines", line = list(width=10), 
#                         showlegend = TRUE,  name='PC1') 
#fig <- fig %>% add_trace(x=~x[2,] , y=~y[2,], z=~z[2,], type="scatter3d", mode="lines", line = list(width=10), showlegend = FALSE) 
#fig <- fig %>% add_trace(x=~x[3,] , y=~y[3,], z=~z[3,], type="scatter3d", mode="lines", line = list(width=10), showlegend = FALSE) 

fig <- fig %>% layout( xaxis = list(title='x'), yaxis = list(title = 'Y Axis Title'), zaxis=list(title = 'z'))
fig

# ----

# Centroid

x <- c(PCA$center[1], PCA$center[1] + PCA$rotation[2,1]*PCA$sdev[1])
y <- c(PCA$center[2], PCA$center[2] + PCA$rotation[2,2]*PCA$sdev[1])
z <- c(PCA$center[3], PCA$center[3] + PCA$rotation[2,3]*PCA$sdev[1])
data_PCA <- data.frame(x, y, z)

# -----

library(plotly)


data <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/gapminderDataFiveYear.csv")


data_2007 <- data[which(data$year == 2007),]

data_2007 <- data_2007[order(data_2007$continent, data_2007$country),]

data_2007$size <- data_2007$pop

colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')


fig <- plot_ly(data_2007, x = ~gdpPercap, y = ~lifeExp, z = ~pop, color = ~continent, size = ~size, colors = colors,
               
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 10),
               
               text = ~paste('Country:', country, '<br>Life Expectancy:', lifeExp, '<br>GDP:', gdpPercap,
                             
                             '<br>Pop.:', pop))

fig <- fig %>% layout(title = 'Life Expectancy v. Per Capita GDP, 2007',
                      
                      scene = list(xaxis = list(title = 'GDP per capita (2000 dollars)',
                                                
                                                gridcolor = 'rgb(255, 255, 255)',
                                                
                                                range = c(2.003297660701705, 5.191505530708712),
                                                
                                                type = 'log',
                                                
                                                zerolinewidth = 1,
                                                
                                                ticklen = 5,
                                                
                                                gridwidth = 2),
                                   
                                   yaxis = list(title = 'Life Expectancy (years)',
                                                
                                                gridcolor = 'rgb(255, 255, 255)',
                                                
                                                range = c(36.12621671352166, 91.72921793264332),
                                                
                                                zerolinewidth = 1,
                                                
                                                ticklen = 5,
                                                
                                                gridwith = 2),
                                   
                                   zaxis = list(title = 'Population',
                                                
                                                gridcolor = 'rgb(255, 255, 255)',
                                                
                                                type = 'log',
                                                
                                                zerolinewidth = 1,
                                                
                                                ticklen = 5,
                                                
                                                gridwith = 2)),
                      
                      paper_bgcolor = 'rgb(243, 243, 243)',
                      
                      plot_bgcolor = 'rgb(243, 243, 243)')


fig



# -----
library(plotly)


mtcars$am[which(mtcars$am == 0)] <- 'Automatic'

mtcars$am[which(mtcars$am == 1)] <- 'Manual'

mtcars$am <- as.factor(mtcars$am)


fig <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E'))

fig <- fig %>% add_markers()

fig <- fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
                                   
                                   yaxis = list(title = 'Gross horsepower'),
                                   
                                   zaxis = list(title = '1/4 mile time')))


fig
