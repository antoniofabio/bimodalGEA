cp DESCRIPTION bimodalGEA
echo Date: `date +%Y-%m-%d` >> bimodalGEA/DESCRIPTION
R CMD build bimodalGEA
