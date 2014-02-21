R --vanilla -f "E:\\sqh\\programs\\r\\papers\\2014_metavalid\\buildpackage.R"
E:
cd "E:\\sqh\\programs"
R CMD check metavalid
R CMD build metavalid
R CMD INSTALL --build metavalid
Rscript -e "install.packages('metavalid_0.99.0.zip', repos = NULL)"

