/***************************************************************************************************************************************
Para computarmos os níveis de produtividade iniciais utilizamos  a equação (19) da seção 5.1 de Costinot at al (2012), medindo assim
a produtividade "revelada" 

Primeiro definindo as variáveis e parâmetros necessários: 

x_ijk = exportação total do país i para  o país j na industria k
delta_ij = exportador-importador (FEs)  
delta_jk = importador-industria (FEs)
delta_ik = exportador-industria (FEs)

Equação (19) pg. 21: 
ln(x_ijk) = delta_ij + delta_jk + delta_ik + eps_ijk

****************************************************************************************************************************************/
cd "~\Área de Trabalho\trade_projeto_1\"

*Base WIOT 2011- organizada 
use wtoi2011.dta, clear

/*
Alterações necesssárias na base WIOT 2011 depois de utilizar o do file disponível para organizá-la no formato stata  
*/

rename out_country expor
rename in_country impor
rename out_product indus
rename Xio x_ijk

*drop Row e TOT
drop if expor=="RoW" | expor == "TOT"
drop if impor == "RoW"


* Construindo x_ijk: gastos do país i importando do país j no setor j, agrega a base para importadores e depois constroi o seu log:

collapse (sum) x_ijk, by (expor impor indus)

*tratando os negativos
replace x_ijk =0 if x_ijk < 0



gen lnx_ijk = ln(x_ijk)

sort expor indus impor

foreach varr in "expor" "impor" "indus" {
rename `varr' `varr'_label
encode `varr'_label, gen(`varr')
}
label drop expor impor indus

export delimited expor expor_label impor impor_label indus indus_label using labelbook.csv, nolabel replace

*Construindo as dummies para cada nível de exportador(i), importador(j) e industria(k):
tab expor, generate(it)
tab impor, generate(jt)
tab indus, generate(kt)


**Construindo agora as dummies necessárias para a eq. 19:

set matsize 11000

*delta_ij = exportador-importador (FEs)  
foreach x of varlist it*{
foreach y of varlist jt*{
generate delta1`x'`y' = `x'*`y'
}
}
*delta_jk = importador-industria (FEs)
foreach x of varlist it*{
foreach y of varlist kt*{
generate delta2`x'`y' = `x'*`y'
}
}
*delta_ik = exportador-industria (FEs)
foreach x of varlist kt*{
foreach y of varlist jt*{
generate delta3`x'`y' = `x'*`y'
}
}		


eststo: qui reg lnx_ijk delta1* delta2* delta3*, nocons
esttab using temp4.csv, replace not nostar b(a4)



*NORMALIZAÇÃO 
insheet using temp4.csv, clear


		rename v1 name
		rename v2 value
		drop if name=="" | name=="N"
		replace value = "" if value=="0"
		destring value, replace
		gen name1 = subinstr(name,"delta1", "",.)
		replace name1 = subinstr(name1, "delta2", "", .)
		replace name1 = subinstr(name1, "delta3", "", .)
		replace name1 = subinstr(name1, "o.", "", .)
		split name1, gen(name1_) parse("t")
		rename name1_1 idtype_1
		rename name1_3 id_2
		gen idtype_2 = substr(name1_2,-1,1)
		gen id_1 = subinstr(name1_2, "j", "", .)
		replace id_1 = subinstr(id_1, "k", "", .)
		drop name1_2
		destring id_1 id_2, replace
		drop name name1
		
		keep if idtype_1=="i" & idtype_2=="k"
		
		rename id_1 i_id
		rename id_2 k_id
			
		drop idtype_1 idtype_2
	
		rename value FE_ik_
		replace FE_ik = 0 if FE_ik ==.

		
		*CAlculando z e normalizando para um em USA e agricultura, usando theta=6.534
		
		gen z = exp(FE_ik_/6.534)
		drop FE_ik_
		rename i_id exporter
		rename k_id industry
		
		gen z_US_agro = z if exporter==40 & industry == 1
		egen z_US_agro_m = mean(z_US_agro)
		drop z_US_agro
		rename z_US_agro_m z_US_agro

		gen z_US = z if exporter==40 
		bys industry:	egen z_US_m = mean(z_US)
		drop z_US
		rename z_US_m z_US
		
		gen z_agro = z if industry == 1
		bys exporter: egen z_agro_m = mean(z_agro)
		drop z_agro
		rename z_agro_m z_agro
				
		gen output = z * z_US_agro/(z_US*z_agro)
		drop z*

		*reshape wide output, i(exporter) j(industry) 			
		*drop exporter
	*		outsheet output* using "zijk2.csv", comma nonames replace
	
	export delimited using zijk2.csv, replace 


/*************Encontrar os parametros da eq(4)h de Costinot at al (2012)******************************
alpha_jk = share de gasto nas variedades da industria k no país j
pi_ijk = share de exportação do país i no país j e industria k
******************************************************************************************************/ 
use wtoi2011.dta, clear

/*
Alterações necesssárias na base WIOT 2011 depois de utilizar o do file disponível para organizá-la no formato stata  
*/

rename out_country expor
rename in_country impor
rename out_product indus
rename Xio x_ijk

*drop Row e TOT
drop if expor=="RoW" | expor == "TOT"
drop if impor == "RoW"


* Construindo x_ijk: gastos do país i importando do país j no setor j, agrega a base para importadores e depois constroi o seu log:

collapse (sum) x_ijk, by (expor impor indus)

*tratando os negativos
replace x_ijk =0 if x_ijk < 0
	
keep expor indus impor x_ijk lnx_ijk

* exportação total por país: 
egen x_j = sum(x_ijk), by(impor)

* exportação total por país em cada setor:
egen x_jk = sum(x_ijk), by(impor indus)

* alpha_jk:
gen alpha_jk = x_jk/x_j 

*pi_ijk:
gen pi_ijk = x_ijk/x_jk 

keep expor indus impor x_ijk alpha_jk pi_ijk
sa alp_pi.dta, replace


*EXPORTANDO A VARIÁVEL LABOR DA wiod_sea já com apenas 2011 nos dados

use wiod_social_2011.dta, clear
keep if Variable == "EMP"
keep if Code == "TOT"
drop Variable Description Code
rename Country expor
rename _2011 L_i
destring L_i, replace
save trabalho.dta

use alp_pi.dta, clear

merge m:1 expor using trabalho.dta

sa vareq4.dta



