################################################################
#         FUNCIONES PARA LOS CALCULOS CON DIAGNOSTICOS         #
################################################################
#
#-------- FUNCIONES PARA MANEJAR DATOS DE COMORBILIDAD ----------------------#
binarias_grupo <- function(dat, dx, grupo, etiqueta = ""){
        # Funcion que coge dat, y para el grupo seleccionado, le "pega" a dat las variables binarias
        # calculadas a partir de dx, distinguiendo por origen
        # Parametros;
        #--dat: base de datos a nivel paciente (esta es la que se rellena con variables)
        #--dx: base de diagnosticos considerados (admite varias filas por paciente)
        #--grupo: el grupo diagnostico seleccionado (un nivel de la variable grupo_propio en dx)
        #--etiqueta: opcional, por si se quiere especificar la fuente (ej.: _hospital)
        # Devuelve:
        # El data frame dat con variables nuevas (binarias y clasificaci?n a nivel paciente)
        # 
        # Filtrar diagnosticos del grupo
        dx <- dx[dx$grupo_propio == grupo, ]
        #
        # Calcular binarias 
        if(nrow(dx) != 0){
                df <- dcast(data = dx %>% mutate(var = 1), id ~ grupo_propio, value.var = "var")
                names(df)[2] <- paste0(grupo, etiqueta)
                nombre1 <- names(df)[2] #lo guardamos
                # Pegarla a dat
                dat <- merge(x = dat, y = df, all.x = TRUE)
                # recodificar NA's a cero:
                dat[,nombre1] <- ifelse(is.na(dat[,nombre1]), 0, dat[,nombre1])
                # Convertirla en binaria: (valores mayores que 1, ser?n 1)
                dat[,nombre1] <- ifelse(dat[,nombre1] > 1, 1, dat[,nombre1])
        } else {# No hay ningun paciente diagnosticado en ese grupo propio
                nombre1 <- paste0(grupo, etiqueta)
                dat[,nombre1] <- rep(0, nrow(dat))
        }
        return(dat)
}
calculo_prevalencias <- function(dat, grupo, comarca = NULL){ #grupo = "anemia", comarca = "COMARCA ARABA"
  # Funcion que calcula prevalencias segun fuentes
  #-- prevalencia del grupo diagnostico (global o por comarca)
  #-- prevalencia de los diagnosticos de AP 
  #-- prevalencia de los diagnosticos hospitalarios
  #-- prevalencia de los diagnosticos de CE 
  #-- porcentaje de pacientes con diagnostico "no redundante" para AP (i.e. solo en AP)
  # Parametros:
  #-- dat: base de datos a nivel paciente con las variables binarias calculadas por binarias_grupo()
  #-- grupo: grupo diagnostico
  #-- comarca: opcional (por defecto calcula los porcentajes globa?es.)
  #
  # Filtramos comarca, si nos la dan
  if(!is.null(comarca)){dat <- dat[dat$comarca == comarca, ]}
  # Crear la variable de grupo.
  var <- grep(names(dat), pattern = paste0(grupo, "_dx"), value = TRUE)
  # Calcular los valores
  #Prevalencia general (sin distinguir fuente)
  #--numero diagnosticados
  n_dx <- sum(dat[, var])
  #--prevalencia Global
  prev_dx <- n_dx / nrow(dat)
  #Prevalencia AP
  #--numero diagnosticados AP
  var1 <- grep(names(dat), pattern = paste0(grupo, "_ap"), value = TRUE)
  n_dx_ap <- sum(dat[, var1])
  #--prevalenia diagnosticados AP
  prev_ap <- n_dx_ap / nrow(dat)
  #Prevalencia hospital
  #--numero diagnosticados hospital
  var2 <- grep(names(dat), pattern = paste0(grupo, "_hosp"), value = TRUE)
  n_dx_hosp <- sum(dat[, var2])
  #--prevalenia diagnosticados hospital
  prev_hosp <- n_dx_hosp / nrow(dat)
  #Prevalencia CE
  #--numero diagnosticados CE
  var3 <- grep(names(dat), pattern = paste0(grupo, "_ce"), value = TRUE)
  n_dx_ce <- sum(dat[, var3])
  #--prevalenia diagnosticados hospital
  prev_ce <- n_dx_ce / nrow(dat)
  #porcentaje de no-redundantes (solo ap), sobre el total 
  varxx <- grep(names(dat), pattern = paste0(grupo, "_ce"), value = TRUE) #ce
  varyy <- grep(names(dat), pattern = paste0(grupo, "_hosp"), value = TRUE) #hosp
  soloAP <- sum(dat[,varxx] == 0 & dat[, varyy] == 0 & dat[, var1] == 1)
  p_soloAP <- soloAP / nrow(dat)
  #
  valores <- c("n_dx" = n_dx, 
               "prev_dx" = round(prev_dx,2), 
               "n_dx_ap" = n_dx_ap, 
               "prev_ap" = round(prev_ap,2), 
               "n_dx_hosp" = n_dx_hosp, 
               "prev_hosp" = round(prev_hosp,2),
               "n_dx_ce" = n_dx_ce, 
               "prev_ce" = round(prev_ce,2),
               "solo_ap" = soloAP,
               "p_solo_ap" = round(p_soloAP,2) )
  return(valores)
}
calculo_Charlson <- function(dat, dx){
  # A la base de datos a nivel paciente (dat), se le pegan las 17 variables binarias
  #correspondientes las componentes del indice de Charlson basandose en la bd de diagnosticos dx
  # Parametros:
  #--dat: base de datos a nivel paciente
  #--dx: base de datos de diagnosticos (id Paciente, codigos CIE)
  # Componente a componente
  dx.acv <- dx %>% filter(cie %in% grep("^430|^431|^432|^433|^434^|^435|^436|^437|^438", x = dx$cie, value = T))
  dat$char_acv <- ifelse(dat$id %in% dx.acv$id, 1, 0)
  rm(dx.acv)
  #
  dx.demencia <- dx %>% filter(cie %in% grep("^290", x = dx$cie, value = T))
  dat$char_demencia <- ifelse(dat$id %in% dx.demencia$id, 1, 0)
  rm(dx.demencia)
  #
  dx.dm_concompli <- dx %>% filter(cie %in% grep("^2504|^2505|^2506", x = dx$cie, value = T))
  dat$char_dm_concompli<- ifelse(dat$id %in% dx.dm_concompli$id, 1, 0)
  rm(dx.dm_concompli )
  #
  dx.dm_sincompli <- dx %>% filter(cie %in% grep("^2500|^2501|^2502|^2503|^2507|\\<250\\>", x = dx$cie, value = T))
  dat$char_dm_sincompli<- ifelse(dat$id %in% dx.dm_sincompli$id, 1, 0)
  rm(dx.dm_sincompli )
  #
  dx.epoc <- dx %>% filter(cie %in% grep("^49|^500|^501|^502|^503|^504|^505|^5064", x = dx$cie, value = T))
  dat$char_epoc<- ifelse(dat$id %in% dx.epoc$id, 1, 0)
  rm(dx.epoc )
  #
  dx.hemiplejia <- dx %>% filter(cie %in% grep("^3441|^342|\\<344\\>", x = dx$cie, value = T))
  dat$char_hemiplejia<- ifelse(dat$id %in% dx.hemiplejia$id, 1, 0)
  rm(dx.hemiplejia )
  #
  dx.iam <- dx %>% filter(cie %in% grep("^410|^412", x = dx$cie, value = T))
  dat$char_iam<- ifelse(dat$id %in% dx.iam$id, 1, 0)
  rm(dx.iam )
  #
  dx.icc <- dx %>% filter(cie %in% grep("^428", x = dx$cie, value = T))
  dat$char_icc<- ifelse(dat$id %in% dx.icc$id, 1, 0)
  rm(dx.icc )
  #
  dx.malignancy <- dx %>% filter(cie %in% grep("^14|^15|^16|^170|^171|^172|^174|^175|^176|^177|^178|^179|^18|
  ^190|^191|^192|^193|^194|^195|^200|^201|^202|^203|^204|^205|
    ^206|^207|^208", x = dx$cie, value = T))
  dat$char_malignancy<- ifelse(dat$id %in% dx.malignancy$id, 1, 0)
  rm(dx.malignancy )
  #
  dx.metastasis <- dx %>% filter(cie %in% grep("^196|^197|^198|^1990|^1991|\\<199\\>", x = dx$cie, value = T))
  dat$char_metastasis<- ifelse(dat$id %in% dx.metastasis$id, 1, 0)
  rm(dx.metastasis )
  # 
  dx.mildliver <- dx %>% filter(cie %in% grep("^5712|^5714|^5715|^5716|\\<571\\>", x = dx$cie, value = T))
  dat$char_mildliver<- ifelse(dat$id %in% dx.mildliver$id, 1, 0)
  rm(dx.mildliver )
  #
  dx.mode_severe_liver <- dx %>% filter(cie %in% grep("^4560|^4561|^4562|^5722|^5723|^5724|^5725|^5726|^5727|^5728|
\\<572\\>", x = dx$cie, value = T))
  dat$char_mode_severe_liver<- ifelse(dat$id %in% dx.mode_severe_liver$id, 1, 0)
  rm(dx.mode_severe_liver )
  #
  dx.peri_vascular <- dx %>% filter(cie %in% grep("^4439|^441|^7854|^V434|\\<443\\>", x = dx$cie, value = T))
  dat$char_peri_vascular<- ifelse(dat$id %in% dx.peri_vascular$id, 1, 0)
  rm(dx.peri_vascular )
  #
  dx.renal <- dx %>% filter(cie %in% grep("^582|^5830|^5831|^5832|^5833|^5834|^5835|^5836|^5837|^585|
^586|^588|\\<583\\>", x = dx$cie, value = T))
  dat$char_renal<- ifelse(dat$id %in% dx.renal$id, 1, 0)
  rm(dx.renal )
  #
  dx.reuma <- dx %>% filter(cie %in% grep("^7100|^7101|^7104|^7140|^7141|^7142|^71481|^725|\\<710\\>|
\\<714\\>", x = dx$cie, value = T))
  dat$char_reuma<- ifelse(dat$id %in% dx.reuma$id, 1, 0)
  rm(dx.reuma )
  #
  dx.ulcera <- dx %>% filter(cie %in% grep("^531|^532|^533|^534", x = dx$cie, value = T))
  dat$char_ulcera<- ifelse(dat$id %in% dx.ulcera$id, 1, 0)
  rm(dx.ulcera )
  #
  dx.vih <- dx %>% filter(cie %in% grep("^042|^043|^044", x = dx$cie, value = T))
  dat$char_vih<- ifelse(dat$id %in% dx.vih$id, 1, 0)
  rm(dx.vih )
  #
  # # Charlson Global
  dat$charlson <- dat$char_acv +
          dat$char_demencia +
          dat$char_dm_sincompli +
          dat$char_epoc +
          dat$char_iam +
          dat$char_icc +
          dat$char_mildliver +
          dat$char_peri_vascular +
          dat$char_reuma +
          dat$char_ulcera +
          dat$char_dm_concompli * 2 +
          dat$char_hemiplejia * 2 +
          dat$char_malignancy * 2 +
          dat$char_renal * 2 +
          dat$char_mode_severe_liver * 3 +
          dat$char_metastasis * 6 +
          dat$char_vih * 6
  # 
  return(dat)
}
comorbilidades_grupo <- function(dat, grupo = "dm", etiqueta = "", vartot = "", what = "porcentaje"){
  # Calcula el porcentaje de pacientes que solo tienen la comorbilidad del grupo (what = "porcentaje)
  #o bien, el numero medio de otras patologias adicionales asociadas al grupo (what = "media")
  # Parametros:
  #-- dat: base de datos a nivel paciente; ha de tener total diagnosticos distintos y total diagnosticos por grupo
  #-- grupo: patologia
  #-- etiqueta: Por si se especifica la fuente (ejemplo:_hospi)
  #-- vartot: Nombre nombre de la variable con el total
  #-- what: porcentaje: % de pacientes; media: numero medio de patologias adicionales; se: error estandar de la media
  #
  # Filtrar pacientes con el grupo
  var <- grep(paste0(grupo, etiqueta), names(dat), value =T)
  dat.grupo <- dat[dat[,var] > 0,]
  if(what == "porcentaje"){
    n_sinotra <- sum(dat.grupo[,vartot] == 1)
    pct_sinotra <- n_sinotra / nrow(dat.grupo)
    return(pct_sinotra)
  }
  if(what == "media"){
    #restamos 1 al total de patologias
    dat.grupo[,"tot_1"] <- dat.grupo[, vartot] -1 
    media_patol <- mean(dat.grupo$tot_1)
    return(media_patol)
  }
  if(what == "se"){
    #restamos 1 al total de patologias
    dat.grupo[,"tot_1"] <- dat.grupo[, vartot] -1 
    sd_patol <- sd(dat.grupo$tot_1)/sqrt(length(dat.grupo$tot_1))
    return(sd_patol)
  }
}
prevalencias_comarca <- function(dat, grupos, comarca = "URIBE", etiqueta = "" ){
        # Devuelve prevalencias de la comarca en un data frame, para una lista de patologias
        # Parametros:
        #-- dat, fichero de datos con variables binarias
        #-- comarca, comarca de interes
        #-- gruposCCS: vector de patologias de las cuales queremos calcular la prevalencia en la comarca
        #-- etiqueta: para poder pasar nombre completo de las variables binarias de patologias (ej. "_hospi")
        df <- as.data.frame(grupos)
        #filtrar comarca
        dat.comarca <- dat[dat$comarca == comarca,]
        preval <- sapply(X = 1:length(grupos), FUN = function(x){
                round(sum(dat.comarca[,paste0(grupos[x], etiqueta)]) / nrow(dat.comarca) * 100, 1)})
        df$comarca <- comarca
        df$preval <- preval
        return(df)
}
#-------- FUNCIONES IDENTIFICACION HOSPITALIZACIONES POTENCIALMENTE EVITABLES ----------------------#
identificar_hpe <- function(dat, label_var = "ing_pe_icc_30d", df_ing, df_ing_sec, df_proc, cie_dg, cie_dg_sec_excluir, cie_proc_excluir){
        # Coge el df a nivel paciente, dat, y le anyade una columna nueva indicando si el paciente ha tenido o no una hospitalizacion PE
        # en el periodo indicado, y en la patologia indicada (a traves de los cies)
        # Parametros:
        #-- dat: df a nivel paciente (dat0)
        #-- label_var: etiqueta para nombrar a la variable nueva (preferible indicar patologia y periodo)
        #-- df_ing: data frame con ingresos de pacientes de dat, dg principales, en el periodo indicado por label_var
        #-- df_ing_sec: data frame con ingresos de pacientes de dat, dg secundarios, en el periodo indicado por label_var
        #-- df_proc: data frame con procedimientos a pacientes de dat, sin filtrar periodo: se identificaran los correspondientes al ingreso indice
        #-- cie_dg: codigos CIE en la identificacion de la patologia indicada, dg principales (introducir manualmente)
        #-- cie_dg_sec_excluir: codigos CIE de dg secundarios para excluir hpe en la patologia indicada (introducir manualmente)
        #-- cie_proc_excluir: codigos CIE de procedimientos para excluir hpe en la patologia indicada (introducir manualmente)
        #
        # dat = dat0
        # label_var = "ing_pe_icc_30d"
        # df_ing = ing_post_30d
        # df_ing_sec = ing_sec_post_30d
        # df_proc = proc
        # cie_dg = ingresos_ic
        # cie_dg_sec_excluir = excluir_ic
        # cie_proc_excluir = excluir_ic_proc
        #
        # 1-Identificar ingresos patologia
        df_ing$ingreso_pat <- 0
        df_ing[grep(pattern = paste(cie_dg, collapse = "|"), x = df_ing$cie),"ingreso_pat"] <- 1
        # 2-Identificar ingresos secundarios que generan exclusion, en base a cies
        df_ing_sec$excluir_hpe <- 0
        df_ing_sec[grep(pattern = paste(cie_dg_sec_excluir, collapse = "|"), x = df_ing_sec$cie), ]$excluir_hpe <- 1
        # 3-Identificar procedimientos secundarios que generan exclusion, en base a cies
        df_proc$excluir_hpe <- 0
        df_proc[grep(pattern = paste(cie_proc_excluir, collapse = "|"), x = df_proc$cie), ]$excluir_hpe <- 1
        #
        # 4-Identificar ingresos PE:
        #
        # ----- PARTE A) Excluir por diagnosticos secundarios
        df <- df_ing_sec %>% group_by(id, ffin) %>% summarise(ind = sum(excluir_hpe)) 
        df$ind <- ifelse(df$ind > 0, 1, 0)
        #nos quedamos con filas con indicador de excluir
        df <- df %>% filter(ind == 1)
        # 2.- en segundo lugar, llevamos este indicador a ingresos principales
        df2 <- merge(x = df_ing[, c("id", "ffin","ingreso_pat")], y = df, by = c("id", "ffin"), all.x = TRUE)
        df2[is.na(df2$ind),"ind" ] <- 0
        #
        # ----- PARTE B) Excluir por Procedimientos
        # i.- A cada procedimiento con criterio de exclusion debemos buscarle su fecha del ingreso principal, para cada paciente
        # Para aligerar el calculo nos quedamos con los que ya provocan exclusion
        df_proc <- df_proc %>% filter(excluir_hpe == 1)
        # Para cada paciente y fecha de realizacion del procedimiento, buscamos la fecha de fin del ingreso correspondiente, con ello creamos una columna nueva en dfproc_post.
        fecvalues <- df_proc$frealizacion
        for(x in 1:nrow(df_proc)){#x=1
                idpac <- df_proc[x,]$id 
                fproc <- df_proc[x,]$frealizacion
                fecvalues[x] <- match_proc_ingreso(id = idpac, frealizacion = fproc, df_ing = df_ing)
        }
        df_proc$ffin <- fecvalues
        # De manera que, en este data frame, tenemos los ingresos (diagnosticos ppales) a los que se les aplicaria la exclusion segun el tipo de procedimiento
        # 1.- para cada paciente y fecha, sumamos el numero de filas en procedimientos a excluir
        df <- df_proc %>% group_by(id, ffin) %>% summarise(ind = sum(excluir_hpe)) 
        df$ind <- ifelse(df$ind > 0, 1, 0)
        #nos quedamos con filas con indicador de excluir
        df <- df %>% filter(ind == 1)
        # 2.- en segundo lugar, llevamos este indicador a ingresos principales
        df2.1 <- merge(x = df_ing[, c("id", "ffin","ingreso_pat")], y = df, by = c("id", "ffin"), all.x = TRUE)
        df2.1[is.na(df2.1$ind),"ind" ] <- 0
        table(df2.1$ind)
        table(df2.1$ingreso_pat, df2.1$ind)
        # Quedarnos con los que cumplirian el criterio de "ingreso ic potencialmente evitable"
        # Para ello, tenemos que mezclar ingresos y procedimientos en un solo data frame
        df22 <- merge(x = df2, y = df2.1, by = c("id", "ffin", "ingreso_pat"), all.x = TRUE)
        df22$ind <- ifelse(df22$ind.x == 1 | df22$ind.y == 1, 1, 0)
        # ----- PARTE C) Calculo final
        # Quedarnos con los que cumplirian el criterio de "ingreso ic potencialmente evitable" segun los dos criterios
        df2final <- df22 %>% filter(ingreso_pat == 1 & ind == 0) #ic, y no excluye ningun criterio.
        # df2 contiene los ingresos posteriores, motivo ppal icc que SON evitables (i.e. NO cumplen ninguna exclusion)
        names(df2final)
        # Indicador de ingreso potencialmente evitable
        # Calcular el numero por id
        df3 <- df2final %>% group_by(id) %>% summarise(num_ing_post_pat_pe = n())
        # Pasar a dat
        dat <- merge(x = dat, y = df3, all.x = TRUE)
        names(dat)[which(names(dat) == "num_ing_post_pat_pe")] <- label_var
        dat[is.na(dat[,label_var]), label_var] <- 0
        #convertirla en binaria
        dat[, label_var] <- ifelse(dat[,label_var] > 0, 1, 0)
        return(dat)
}
match_proc_ingreso <- function(id, frealizacion, df_ing){#id = 517921379 ; frealizacion = "2011-10-24"
        #funcion auxiliar que devuelve el ingreso que corresponde al procedimiento
        idpac <- id; fpac <- frealizacion
        #filtrar ingresos principales del paciente, con fecha inicio posterior al proc
        dfing <- df_ing %>% filter(id == idpac, fpac >= finicio, fpac <= ffin)
        #deberia tener una unica fila.
        ffin <- dfing$ffin
        #NOTA. Hay casos en que un mismo procedimiento encaja entre las fechas de dos ingresos (finicio, ffin); 
        #como no tenemos manera de saber a cual corresponde, cogemos siempre el segundo
        if(length(ffin) > 1){
                ffin <- ffin[2]
        }
        return(ffin)
}
#-------- FUNCIONES AUXILIARES REPRESENTACIONES GRAFICAS---------------#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
ggCaterpillar <- function(re, QQ=TRUE, likeDotplot=TRUE) {
        require(ggplot2)
        f <- function(x) {
                pv   <- attr(x, "postVar")
                cols <- 1:(dim(pv)[1])
                se   <- unlist(lapply(cols, function(i) sqrt(pv[i, i, ])))
                ord  <- unlist(lapply(x, order)) + rep((0:(ncol(x) - 1)) * nrow(x), each=nrow(x))
                pDf  <- data.frame(y=unlist(x)[ord],
                                   ci=1.96*se[ord],
                                   nQQ=rep(qnorm(ppoints(nrow(x))), ncol(x)),
                                   ID=factor(rep(rownames(x), ncol(x))[ord], levels=rownames(x)[ord]),
                                   ind=gl(ncol(x), nrow(x), labels=names(x)))
                
                if(QQ) {  ## normal QQ-plot
                        p <- ggplot(pDf, aes(nQQ, y))
                        p <- p + facet_wrap(~ ind, scales="free")
                        p <- p + xlab("Standard normal quantiles") + ylab("Random effect quantiles")
                } else {  ## caterpillar dotplot
                        p <- ggplot(pDf, aes(ID, y)) + coord_flip()
                        if(likeDotplot) {  ## imitate dotplot() -> same scales for random effects
                                p <- p + facet_wrap(~ ind)
                        } else {           ## different scales for random effects
                                p <- p + facet_grid(ind ~ ., scales="free_y")
                        }
                        p <- p + xlab("Levels") + ylab("Random effects")
                }
                
                p <- p + theme(legend.position="none")
                p <- p + geom_hline(yintercept=0)
                p <- p + geom_errorbar(aes(ymin=y-ci, ymax=y+ci), width=0, colour="black")
                p <- p + geom_point(aes(size=1.2), colour="blue")
                return(p)
        }
        
        lapply(re, f)
}
# --------- FUNCIONES OBSOLETAS ---------------- #
# #-------- FUNCIONES PARA CREAR ESTRUCTURA PARA  MODELO MIXTO AP -----------#
rellenar_var <- function(grupo, df){#grupo_propio = "dm"
  # Al data frame df (estructura id * grupo_propio y variable var), rellena var con un 1 si
  # el paciente, en el grupo "grupo", ha sido diagnosticado en AP.
  # Si el paciente NO ha sido diagnosticado en AP en el grupo, entonces mantiene la variable original
  # - que indica si ha sido diagnosticado por alguna de las fuentes (=0) o no (=NA)
  # Par?metros:
  ## grupo diagnostico
  ## df: data frame con estructura indicada.
  #seleccionar filas del grupo en pdf
  var_ap <- paste0(grupo, "_ap")
  #en dat0, buscamos si est? diagnosticado en AP
  dat0 <- dat0[dat0[,var_ap] == 1, ]
  df <- df %>% mutate(ind = ifelse(id %in% dat0$id & grupo_propio == grupo, 1, 0))
  df[df$ind == 1, "var"] <- 1
  #
  return(df)
}
# #-------- FUNCIONES PARA REPRESENTACI?N GR?FICA ---------------------------#
# dibujar_grupo_comarcas <- function(dat, grupo = "dm", fuente = "ap", redundantes = "no"){
#   # Funci?n que representa prevalencias por cada grupo diagnostico, por comarcas
#   # Par?metros:
#   #-- dat: base de datos (1 fila por paciente), con las variables binarias: _ap, _ce, _dx
#   #donde: _dx vale 1 si el paciente ha sido diagnosticado en cualquiera de las fuentes; si no, vale 0
#   #donde: _ap vale 1 si el paciente ha sido diagnosticado en AP, anteriormente al ingreso, y no se ha cerrado antes; si no, vale 0
#   #donde: _ce vale 1 si el paciente ha sido diagnosticado en AP, anteriormente al ingreso, y no se ha cerrado antes; si no, vale 0
#   #-- grupo: grupo diagnostico (patolog?a). Por defecto: "dm"
#   #-- fuente: ap o ce (subregistros posibles). Por defecto: "ap"
#   #
#   ## Primer dato a representar: prevalencia en general
#   #La variable ya nos viene calculada en dat.
#   nombre <- paste0(grupo, "_dx")
#   tit <- paste0("Grupo: ", grupo) #titulo del grafico.
#   ## Segundo dato a representar: prevalencia AP
#   # Calculamos df auxiliar
#   col_name <- as.name(nombre)
#   col_name1 <- as.name(paste0(grupo, "_ap"))
#   if(fuente == "ce"){ 
#     col_name1 <- as.name(paste0(grupo, "_ce"))
#   }
#   # Agregar por comarcas, si no nos interesan los redundantes.
#   df <- dat %>% group_by(comarca) %>% summarise(n = n(),
#                                                 s1 = sum(col_name), 
#                                                 pct1 = s1 / n,
#                                                 s2 = sum(col_name1),
#                                                 pct2 = s2 / n)
#   # Si incluye redundantes, los calculamos para el grupo seleccionado
#   if(redundantes == "si"){
#     #calculamos la variable que indica si el paciente se ha diagnosticado en AP, pero NO en hospital
#     var <- paste0(grupo, "_ap_no_hospital")
#     dat[, var] <- ifelse(dat[,paste0(grupo, "_ap")] == 1 & dat[,paste0(grupo, "_hosp")] == 0 , 1, 0)
#     col_name2 <- as.name(var)
#     #agregar por comarcas (para representar)
#     df <- dat %>% group_by(comarca) %>% summarise(n = n(),
#                                                   s1 = sum(col_name), 
#                                                   pct1 = s1 / n,
#                                                   s2 = sum(col_name1),
#                                                   pct2 = s2 / n,
#                                                   s3 = sum(col_name2),
#                                                   pct3 = s3 / n)
#   }
#   # Representar, seg?n AP o CE.
#   if( fuente == "ap" ){#valores de AP en rojo
#     p <- ggplot(data = df, aes(x = comarca, y = pct1)) + 
#       geom_point(stat = "identity", col = "blue", size = 3) + 
#       geom_line(stat = "identity", col = "blue", aes(group = 1),linetype = 2)+
#       coord_flip() + 
#       geom_text(aes(label = round(pct1*100,1)), size = 3, col = "blue", hjust = -0.5) +
#       # ggtitle(tit) +
#       # ylab("% pacientes diagnosticados\n(sin distinguir origen)") +
#       scale_y_continuous(limits = c(0,1)) +
#       #   geom_hline(aes(yintercept=mean(df$pct1)),linetype = "longdash", col ="blue") +
#       geom_point(aes(x = comarca, y = pct2),stat = "identity", col = "#CC0000", size = 3) +
#       geom_line(stat = "identity", col = "#CC0000", aes(x = comarca, y = pct2, group = 1), linetype = 2) +
#       # geom_text(aes(label = round(pct2*100,1)), size = 3, col = "orange", hjust = 2) +
#       ylab("") +
#       ggtitle(tit) 
#     # +
#     # geom_hline(aes(yintercept=mean(df$pct2)),linetype = "longdash", col ="#D55E00")
#     if( redundantes == "si" ){
#       p <- p + geom_line(aes(x = comarca, y = pct3), col = "#CC0000", group = 1, linetype = 15)
#     }
#   }
#   # 
#   if( fuente == "ce" ){ # Valores de CE en verde.
#     p <- ggplot(data = df, aes(x = comarca, y = pct1)) + 
#       geom_point(stat = "identity", col = "blue", size = 3) + 
#       geom_line(stat = "identity", col = "blue", aes(group = 1), linetype = 2)+
#       coord_flip() + 
#       geom_text(aes(label = round(pct1*100,1)), size = 3, col = "blue", hjust = -0.5) +
#       # ggtitle(tit) +
#       # ylab("% pacientes diagnosticados\n(sin distinguir origen)") +
#       scale_y_continuous(limits = c(0,1)) +
#       #   geom_hline(aes(yintercept=mean(df$pct1)),linetype = "longdash", col ="blue") +
#       geom_point(aes(x = comarca, y = pct2),stat = "identity", col = "green", size = 3) +
#       geom_line(stat = "identity", col = "green", aes(x = comarca, y = pct2, group = 1), linetype = 2) +
#       # geom_text(aes(label = round(pct2*100,1)), size = 3, col = "orange", hjust = 2) +
#       ylab("") +
#       ggtitle(tit) 
#     #
#     if( redundantes == "si" ){
#       p <- p + geom_line(aes(x = comarca, y = pct3), col = "#CC0000", group = 1, linetype = 15)
#     }
#   }
#   #
#   # A?adir capa de redundantes, si se ha elegido la opci?n
#   return(p)
# }
# ## re = object of class ranef.mer
#-------- FUNCIONES PARA CALCULAR OUTCOMES ARTICULO 1 E INFORME -----------#
ingresos_servicio <- function(dat, ing, servicio){#servicio = "RESPIRATORIO" ; ing = ing_1a
        # Funcion que coge el df de ingresos ing, y para el servicio seleccionado, le "pega" a dat las variables binarias
        # calculadas a partir de ing
        # Parametros;
        #--dat: base de datos a nivel paciente (esta es la que se rellena con variables)
        #--ing: base de datos con ingresos hospitalarios (con id y servicio) - puede ser ing_prev, ing_1a
        #--servicio: el servicio seleccionado
        # Devuelve:
        # El data frame dat con variables nuevas (binarias)
        # 
        # Filtrar diagnosticos del grupo
        ing <- ing[ing$servicio == servicio, ]
        #
        if(nrow(ing)!=0){
                df <- ing %>% group_by(id) %>% summarise(n_ing = n())
        }
        # cambiar nombre de la variable
        varname <- paste0("n_ing_", tolower(servicio))
        names(df)[2] <- varname
        # llevar a dat
        dat <- merge(x = dat, y = df, all.x = TRUE)
        #recodificar a 0
        dat[is.na(dat[,varname]),varname] <- 0
        return(dat)
}