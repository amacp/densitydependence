library(Rcpp)
paste(getwd(),"/Cppfiles", sep="")
# Make R package
Rcpp.package.skeleton(
    name = "DDS", 
    code_files = c("Rfiles/ShinyApp/server.r", "Rfiles/ShinyApp/ui.r"),
    cpp_files = c(
    "Cppfiles/main.cpp",
    "Cppfiles/random.h", 
    "Cppfiles/random.cpp", 
    "Cppfiles/utils.h", 
    "Cppfiles/utils.cpp",
    "Cppfiles/Makevars"), 
    example_code = FALSE,
    author = "F.J.H. de Haas &  A.M. Pearson",
    email = "dehaas@zoology.ubc.ca")
    
#Add to the DESCRIPTION file:
#Imports: Rcpp (>= 0.12.17), BH, RcppProgress
#LinkingTo: Rcpp, BH, RcppProgress

