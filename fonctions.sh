#=========
#Fonctions
#=========

#Affiche un text avec horodatage
affiche_txt()
{
    horod=$(date '+%D-%Hh%Mm%Ss')
    echo -e $horod"\t"$*
}

#Fonction pour tester les erreurs
test_errors()
{
    errorcode=$1
    errmsg=$2
    progexit=$3
    if [ $errorcode -gt 0 ]
    then
	affiche_txt "ERROR $errorcode: "$errmsg
	if [ $progexit -eq 1 ]
	then
	    exit $errorcode
	fi
    fi
}

#fonction pour tester si un fichier existe et si il n'est pas vide
test_path()
{
    if [ -f "$1" ]
    then
	if [ -s $1 ]
	then
	    return 1 #le fichier existe et il n'est pas vide
	else
	    return 2 #le fichier existe mais il est vide
	fi
    else
	return 0 #le fichier n'existe pas
    fi
}
