# In caso...pulisco -------------------------------------------------------

# rm(list = ls())

# Creo un array di tabelle ------------------------------------------------

# Inizializzo
mie.table = array(NA, c(3,4,20))
# mie.table
# mie.table[,,1]

# Riempio
for(i in 1:dim(mie.table)[3]) mie.table[,,i] = matrix(sample(1:20,3*4,rep = T), 3, 4)
# mie.table

# Ripeto la prima 3 volte e la dodicesima 2 volte
mie.table[,,5]  = mie.table[,,1]
mie.table[,,10] = mie.table[,,1]
mie.table[,,15] = mie.table[,,1]

mie.table[,,7]  = mie.table[,,12]
mie.table[,,14] = mie.table[,,12]


# Confronto ---------------------------------------------------------------

compaTab = function(mie.table){
  eqv.tab = list()                     # Inizializzo il listone contenente gli indici delle tabelle equivalenti...
  eqv.idx = 0                          # ...e l'indice che le scorre...
  act.tab = 1:dim(mie.table)[3]        # Inizializzo indice tabelle attive
  
  while ( length(act.tab) > 0 ){                            # Finché ho tabelle attive...
    ref.tab = act.tab[1]                                    # ...seleziono la prima attiva come tabella di riferimento
    eqv.idx = eqv.idx + 1                                   # ...aggiungo un elemento alla lista delle tabelle equivalenti...  
    eqv.tab[[eqv.idx]] = ref.tab                            # ...lo inizializzo...
    names(eqv.tab)[eqv.idx] <- paste("Table", ref.tab)      # ...lo battezzo...
    act.tab = act.tab[-1]                                   # ...lo rimuovo dalla lista della tabelle attive
    act.idx = act.tab                                       # Salva gli indici per il ciclo <for> successivo... 
    for (j in act.idx){                                     # Confronto la prima ancora attiva con tutte le successive ancora attive...
      if ( all(mie.table[,,ref.tab] == mie.table[,,j]) ){   # ...se uguale...
        eqv.tab[[eqv.idx]] = c(eqv.tab[[eqv.idx]], j)       # ...aggiungo l'indice attuale al mucchio...
        act.tab = setdiff(act.tab, j)                       # ...e lo rimuovo da quelle attive...
      }
    }  
  }
  return(eqv.tab)
}



# Do un'occhiata... -------------------------------------------------------

out = compaTab(mie.table)
out

#...e in effetti la prima e la settima stanno assieme alle tabelle giuste...
# Per avere le frequenze...

lapply(out, length)

#...o anche...
as.matrix(lapply(out, length))
