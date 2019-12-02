#' @title plot.
#' @description  plot scatter plot.
#' @param n int.number:such as 1.2
#' @importFrom graphics plot
#' @export
#' @examples
#' p<-test_plot(5)
#' p
test_plot<-function(n){
  plot(1:n)
}


#' @title test_ggplot
#' @param n input a positive interger number
#' @import ggplot2
#' @export
#' @examples
#' test_ggplot(12)
test_ggplot<-function(n){
  #library(ggplot2)
  ggplot()+geom_point(aes(x=1:n,y=1:n,col="orange"))+
    labs(x="x",y="y",title="test for ggplot",caption = "@qliu")+
    theme(legend.position = "none")
}


