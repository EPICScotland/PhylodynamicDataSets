# functions to classify bird names to latin species names and countries to continents
# S. J. Lycett
# 14 Dec 2011 - extract from from matrixTranslate.R
# 3 April 2012 - added great cormorant to pel, and Bulgaria to Europe
# 23 May 2012 - added more birds
# 24 May 2012 - added more countries
# 13 July 2012 - changed egret from ans to cic
# 5 Dec 2012 - more countries
# 25 Mar 2013 - more birds and Iceland in Europe
# 5 April 2013 - more birds, and corrected pel and cic orders (egret is a pel, not cic)
# 9 April 2013 - added countryClass - change place to country
# 25 April 2013 - getWildDomestic
# 26 April 2013 - problem with grep and owl
# 1  May   2013 - domestic_green-winged_teal => Domestic
# 7  May   2013 - Wenzhou & Zhangzhou in China
# 8  May   2013 - Cuckoo & abbreviations (mostly for HKU)
# 9  May   2013 - CU = chukar
# 5  Aug   2013 - african stonechat
# 2  Dec   2013 - finch (pas)
# 12 Feb   2014 - dove (col)
# 26 Feb   2014 - bantam = chicken = domestic gal, silky fowl = silkie chicken = domestic gal

######################################################################
# function to classify birds in H5N1 - to bird order

# Bird Orders

# Acc	Accipitriformes	http://en.wikipedia.org/wiki/Accipitriformes	Eagles, Hawks (falcons are separate)
# Ans	Anseriformes	http://tolweb.org/Anseriformes	Ducks, Geese, Swans
# Cha	Charadriiformes	http://tolweb.org/Charadriiformes	Shorebirds
# Cic	Ciconiiformes	http://tolweb.org/Ciconiiformes/26332	Storks, Herons
# Cor Coraciiformes	http://en.wikipedia.org/wiki/Coraciiformes - like pas.
# Col Columbiformes	http://en.wikipedia.org/wiki/Columbiformes, doves and pigeons
# Cuc Cuculiformes	http://en.wikipedia.org/wiki/Cuculiformes	Cuckoos
# Gal	Galliformes	http://tolweb.org/Galliformes	Fowl, Quail
# Gru	Gruiformes	http://en.wikipedia.org/wiki/Gruiformes	Coot, Moorhen
# Fal	Falconiformes	http://tolweb.org/Falconiformes/26379	Falcons
# Pas	Passeriformes	http://tolweb.org/Passeriformes	Perching birds (e.g Sparrows)
# Pel Pelecaniforms	http://en.wikipedia.org/wiki/Pelecaniformes Pelican
# Pic Piciformes		http://en.wikipedia.org/wiki/Piciformes - like pas, includes woodpeckers
# Pro	Procellariiformes	http://en.wikipedia.org/wiki/Procellariiformes	Shearwater
# Psi Psittaciformes	http://en.wikipedia.org/wiki/Parrot	Parrot
# Rhe	Rheiformes	http://en.wikipedia.org/wiki/Rhea_(bird)	Rhea (looks like an emu)
# Str	Struthioniformes	http://en.wikipedia.org/wiki/Ostrich	Ostrich, Emu
# Tin Tinamiformes	http://en.wikipedia.org/wiki/Red-winged_Tinamou Tinamou (Brazil)

# birdNames is a vector

# updated 25 Mar 2013
# updated 26 April 2013
# updated 9  May   2013
# updated 12 Feb   2014
birdClass	<- function( birdNames ) {

	birdNames	<- tolower(birdNames)

	# nasty coding of ck for chicken messes up duck !!
	inds		<- which((birdNames=="ck") | (birdNames=="chick"))
	if (length(inds) > 0) {
		birdNames[inds] <- "chicken"
	}	

	inds		<- which(birdNames=="gf")
	if (length(inds) > 0) {
		birdNames[inds] <- "guineafowl"
	}

	inds		<- which(birdNames=="ph")
	if (length(inds) > 0) {
		birdNames[inds] <- "pheasant"
	}

	inds		<- which(birdNames=="qa")
	if (length(inds) > 0) {
		birdNames[inds] <- "quail"
	}


	# crowned gets synoned with crow - so remove all crowneds

	#inds		<- grep("crowned", birdNames)
	#if (length(inds) > 0) {
	#	for (i in 1:length(inds)) {
	#		birdNames[inds[i]] = sub("crowned", "", birdNames[inds[i]])
	#		}
	#}

	# quicker ! not that it matters terribly much, but still
	birdNames <- gsub("crowned", "", birdNames)

	birdOrders	<- array("-", length(birdNames))

	acc	<- c("eagle","hawk","buzzard","bussard","hooded vulture")

	ans	<- c("dk","duck","mallard","teal","goose","gs","swan","cygnus","wigeon",
			"munia","grebe","goldeneye","pochard","goosander",
			"waterfowl","pintail","gadwall","garganey","platyrhynchos","crecca",
			"plathyrhynchos","shoveler","aquatic bird","aquaticbird",
			"widgeon","bufflehead","eider","scoter","redhead","water fowl","watercock","anas acuta",
			"greater scaup","mergus albellus","smew",
			"canvasback","hooded merganser","lesser scaup","northern shoverl","brant")

	cic	<- c("stork","condor")

	cor	<- c("rollers")

	col	<- c("peaceful dove","rock dove", "dove")

	cha	<- c("shorebird","gull","tern","turnstone","dunlin","sandpiper","sanderling",
			"red knot","redknot","knot","plover","red-necked stint","redneckedstint","curlew",
			"rufous-necked stint","black-legged kittiwake","common murre",
			"arenaria interpres","thick-billed murre","arenaria-interpres", "eurasian woodcock")

	cuc	<- c("cuckoo")

	gal	<- c("chicken","pheasant","quail","turkey","sck","guinea","guineafowl","guinea_fowl",
			"chukar","bustard","peacock","peahen","peafowl","partridge",
			"poultry","chukkar","chukka","bantam","silky fowl")

	gru	<- c("coot","moorhen","hooded crane","black-neck crane")

	fal	<- c("falcon","kestrel","peregrine","harrier","owl")

	pas	<- c("sparrow","crow","magpie","pigeon","blackbird","myna",
			"starling","wild bird","wildbird","shrike","robin","softbill",
			"japanese white eye","japanesewhiteeye","common iora","commoniora",
			"fairy bluebird","fairybluebird","rook",
			"babbler","black bulbul","golden mountain thrush",
			"japanese white-eye","silver-eared mesia","brambling","chinese hwamei",
			"african stonechat","finch")

	psi	<- c("parrot","conure","parakeet","psittacine","macaw")

	pro	<- c("shearwater")	
	
	rhe	<- c("rhea")

	str	<- c("ostrich","emu")

	tin	<- c("red-winged tinamou", "redwingedtinamou")

	pel	<- c("pelican","great cormorant","greatcormorant", "cormorant", "ardea cinerea", "grey heron", "egret", "heron")

	pic	<- c("great barbet")

	orders <- list(	acc=acc,ans=ans,cic=cic,cor=cor,col=col,
				cha=cha,cuc=cuc,gru=gru,gal=gal,fal=fal,
				pas=pas,pel=pel,pro=pro,pic=pic,rhe=rhe,str=str,psi=psi,tin=tin)

	for (i in 1:length(orders)) {
		oo	<- unlist(orders[i])
		for (j in 1:length(oo)) {
			if (oo[j] != "owl") {
				inds	<- grep(oo[j], birdNames)
			} else {
				inds  <- which(birdNames==oo[j])
			}
			if (length(inds) > 0) {
				birdOrders[inds] <- attributes(orders[i])$names
			}
		}
	}	

	inds	<- grep("feces", birdNames)
	if (length(inds) > 0) {
		birdOrders[inds] <- "-"
	}

	return( birdOrders )

}

# 2 Sept 2013 - added turkey (not nec. previously), and muscovy_duck as Domestic
# see also http://www.ncbi.nlm.nih.gov/pubmed/20521681 (Qinghai lake and teals)
# During 2007-08 we marked wild ducks at Poyang Lake with satellite transmitters to examine the location and timing of spring migration 
# and identify any spatiotemporal relationship with HPAI H5N1 outbreaks. 
# Species included the 
# Eurasian wigeon (Anas penelope), 
# northern pintail (Anas acuta), 
# common teal (Anas crecca), 
# falcated teal (Anas falcata), 
# Baikal teal (Anas formosa), 
# mallard (Anas platyrhynchos), 
# garganey (Anas querquedula), and 
# Chinese spotbill (Anas poecilohyncha). 
# These wild ducks (excluding the resident mallard and Chinese spotbill ducks) ..

getWildDomestic <- function( bird ) {
	bird2			<- tolower(bird)
	bird2			<- gsub(" ","_",bird2)			# 2 Sept 2013

	wildDomestic	<- array("Wild", length(bird2))
	inds			<- which( 	(bird2=="ck") | 
					(bird2=="chicken") | 
					(bird2=="guinea_fowl") | 
					(bird2=="guineafowl") |
					(bird2=="sck") |
					(bird2=="silky_chicken") |
					(bird2=="silkie_chicken") |
					(bird2=="silkie") |
					(bird2=="silky_fowl") |
					(bird2=="quail") |
					(bird2=="gs") |
					(bird2=="goose") |
					(bird2=="dk") |
					(bird2=="duck") |
					(bird2=="pheasant") |
					(bird2=="ph") |
					(bird2=="partridge") |
					(bird2=="Gf") |
					(bird2=="turkey") |
					(bird2=="village_chicken") |
					(bird2=="muscovy_duck") |
					(bird2=="poultry") |
					(bird2=="bantam") )
	
	wildDomestic[inds] <- "Domestic"
	inds			 <- which(bird2=="mammal")
	wildDomestic[inds] <- "Domestic"
	inds			 <- which(bird2=="human")
	wildDomestic[inds] <- "Domestic"
	inds			 <- which(bird2=="environment")
	wildDomestic[inds] <- "Domestic"
	inds			 <- grep("domestic",bird2)
	wildDomestic[inds] <- "Domestic"

	inds			 <- grep("wild",bird2)
	wildDomestic[inds] <- "Wild"

	return (wildDomestic) 
}

unabbreviatedBird <- function( bird ) {
	bird2		<- toupper(bird)
	abbrev	<- c("CK","CU","DK","GS","PG","PH","SCK","WDK","WWF")
	unabbrev	<- c("chicken","chukar","duck","goose","pigeon","pheasant","chicken","wild duck","wild water fowl")
	minds		<- match(bird2, abbrev)
	jj		<- which(is.finite(minds))
	ubird		<- bird
	ubird[jj]	<- unabbrev[minds[jj]]
	return( ubird )
}

#####################################################################################################
# function to classify countries to continents

# check this before use - it only has limited values
# updated 16 Nov 2011
# updated 31 Jan 2012
# updated 14 Feb 2012
# updated 17 Feb 2012
# updated 24 May 2012
# updated 9  Jul 2012
# updated 5  Dec 2012
# updated 25 Mar 2013
continentClass	<- function( countries ) {
	Europe		<- c(	"UK","France","Austria","Belgium","Bolivia","CzechRepublic","Denmark","Estonia",
					"Germany","Greece","Italy","Luxembourg","Netherlands",
					"Norway","Poland","Turkey","Spain","Serbia",
					"UnitedKingdom","Ireland","Hungary","Sweden","Portugal",
					"Slovenia","Croatia","Herzegovina", "Bulgaria",
					"Switzerland","BosniaandHerzegovina","Slovakia","Finland","Iceland")

	NorthAmerica	<- c("USA","Canada","Mexico")
	CentralAmerica	<- c("ElSalvador","Nicaragua","Haiti","DominicanRepublic",
					"Guatemala","Panama","PuertoRico")

	SouthAmerica	<- c("Chile","Colombia","Argentina","Ecuador","Peru","Brazil","Uruguay")

	Africa		<- c("Ethiopia","Nigeria","Mali","Senegal","Djibouti",
					"Africa","Sudan","Egypt","BurkinaFaso", "Niger", "Kenya",
					"IvoryCoast", "CotedIvoire", "Zambia", "SouthAfrica","Ghana")

	Asia			<- c("India","Malaysia","HongKong","Hong_Kong","China","Singapore","SouthKorea","Thailand","Japan",
					"Cambodia","Taiwan","Myanmar","VietNam","Viet_Nam","Vietnam","Pakistan",
					"Bangladesh","Bhutan", "Indonesia", "Laos","Korea","Malaya",
					"Nepal","Philippines")

	Oceania		<- c("Australia","NewZealand","New_Zealand","Guam","NewCaledonia")

	MiddleEast		<- c("Jordan","Israel","Afghanistan","Arabia",
					"Kuwait","GazaStrip", "Iran", "Iraq", "SaudiArabia", "Lebanon",
					"UnitedArabEmirates","Qatar")

	CentralAsia		<- c(	"Kyrgyzstan","Turkmenistan","Russia",
					"Mongolia","Kazakhstan","Ukraine","Uzbekistan",
					"Georgia","Crimea")

	continents		<- array("-", length(countries))

	conts			<- list( 	Europe=Europe, NorthAmerica=NorthAmerica,
						SouthAmerica=c(CentralAmerica,SouthAmerica),
						Africa=Africa, Oceania=Oceania, Asia=c(Asia,CentralAsia,MiddleEast) )

	for (i in 1:length(conts)) {
		oo 	<- unlist(conts[i])
		
		for (j in 1:length(oo)) {
			inds 	<- which(countries==oo[j])
			if (length(inds) > 0) {
				continents[inds] <- attributes(conts[i])$names
			}	 
		}
	}
	
	return (continents)
}



countryClass	<- function( places ) {

	USA   <- c( "Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut",
			"Delaware","Florida","Georgia","Hawaii","Idaho","Illinois",
			"Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts",
			"Michigan","Minnesota","Mississippi","Missouri",
			"Montana","Nebraska","Nevada","New Hampshire",
			"New Jersey","New Mexico","New York","North Carolina",
			"North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania",
			"Rhode Island","South Carolina","South Dakota","Tennessee",
			"Texas","Utah","Vermont","Virginia",
			"Washington","West Virginia","Wisconsin","Wyoming","Delaware Bay")
	Usa2  <- gsub(" ", "_", USA)
	Usa3	<- gsub(" ", "", USA)
	Usa4  <- c("CA","DE","Memphis","MN","NJ","NY","LA")
	USA   <- sort(unique(c(USA,Usa2,Usa3,Usa4)))

	Canada <- c("Ontario","Quebec","Nova Scotia","New Brunswick",
			"Manitoba","British Columbia","Prince Edward Island","Saskatchewan",
			"Alberta","Newfoundland and Labrador","Newfoundland","Labrador")
	Canada2<- gsub(" ", "_", Canada)
	Canada3<- gsub(" ", "", Canada)
	Canada4<- c("ALB")
	Canada <- sort(unique(c(Canada,Canada2,Canada3,Canada4)))

	Thailand <- c(	"Bangkok","Kanchanaburi","Nakhon Si Thammarat Province","Phra Nakhon Si Ayutthaya","Satun",
				"Amnat Charoen","Khon Kaen","Nan","Phrae","Sing Buri",
				"Ang Thong","Krabi","Narathiwat","Phuket","Sisaket",
				"Bueng Kan","Lampang","Nong Bua Lamphu","Prachinburi","Songkhla",
				"Buriram","Lamphun","Nong Khai","Prachuap Khiri Khan","Sukhothai",
				"Chachoengsao","Loei Province","Nonthaburi","Ranong","Suphan Buri",
				"Chainat","Lopburi Province","Pathum Thani","Ratchaburi","Surat Thani",
				"Chaiyaphum","Mae Hong Son","Pattani","Rayong","Surin",
				"Chanthaburi","Maha Sarakham","Phang Nga","Roi Et","Tak",
				"Chiang Mai","Mukdahan","Phatthalung","Sa Kaeo","Trang",
				"Chiang Rai","Nakhon Nayok","Phayao","Sakon Nakhon","Trat",
				"Chonburi","Nakhon Pathom","Phetchabun","Samut Prakan","Ubon Ratchathani",
				"Chumphon","Nakhon Phanom","Phetchaburi","Samut Sakhon","Udon Thani",
				"Kalasin","Nakhon Ratchasima","Phichit","Samut Songkhram","Uthai Thani",
				"Kamphaeng Phet","Nakhon Sawan","Phitsanulok","Saraburi","Uttaradit",
				"Yala","Yasothon",
				"Phathumthani","Prachinburi_Thailand","NaraThiwat","Nakorn-Patom","Kohn Kaen",
				"Nakhonsawan","Suphanburi" )
	Thailand2 <- gsub(" ", "_", Thailand)
	Thailand3 <- gsub(" ", "", Thailand)
	Thailand  <- sort(unique(c(Thailand,Thailand2,Thailand3)))

	Indonesia <- c("East Java","West Java","East Kalimantan")
	Indonesia2<- gsub(" ","_",Indonesia)
	Indonesia3<- gsub(" ", "", Indonesia)
	Indonesia <- c(Indonesia,Indonesia2,Indonesia3,"Java","Kalimantan")

	China <- c(	"Beijing","Tianjin","Hebei","Shanxi","Mongolia",
			"Liaoning","Jilin","Heilongjiang",
			"Shanghai","Jiangsu","Zhejiang","Anhui","Fujian",
			"Jiangxi","Jiang_Xi","Jiang Xi","Shandong",
			"Henan","Hubei","Hunan","Guangdong","Guangxi",
			"Hainan","Chongqing","Sichuan","Guizhou","Yunnan",
			"Tibet","Xizang","Shaanxi","Gansu",
			"Qinghai","Ningxia","Xinjiang","HongKong","Hong_Kong","HK",
			"Xianggang","Macau","Aomen","Taiwan",
			"Nanchang","Nanjing","Qianzhou","Anyang","Guiyang","Hongze","Huadong","Jiawang",
			"Shantou","Shenzhen","Shijiazhuang","Shuanggou","Taixing","Taizhou","Tongshan",
			"Wuxi","Xianghai","Xiangshui","Xigou","Xuzhou","Yangzhou","Yongcheng","ZhaLong","Dawang","Zibo",
			"jiyuan","Kaifeng","Qixian","Zhuhai","Hangzhou","Wenzhou","Zhangzhou","Guangzhou","Kuming")

	Japan <- c("Kinai","Tokaido","Tosando","Hokurikudo","Sanindo","Sanyodo","Nankaido","Hokkaido",
			"Osaka","Kobe","Yokohama","Akita","Aomori","Chiba","Miyagi","Shiga","Shimane","Tsukuba")

	Korea <- c("Geumgang","Nakdonggang","Seongdong","Shihwa","Hadoree")

	Australia <- c("Victoria","Western Australia","South Australia","Queensland","New South Wales")
	Australia2<- gsub(" ","_",Australia)
	Australia3<- gsub(" ", "", Australia)
	Australia <- sort(unique(c(Australia,Australia2,Australia3)))

	Russia <- c("Astrakhan","Neimonggu","Primorie","Altai","Gurjev","Siberia")

	Germany <- c("Rugen","Berlin","Heinersdorf","Potsdam")

	countries		<- array("-", length(places))

	conts			<- list( USA=USA, Canada=Canada, China=China, Japan=Japan, Thailand=Thailand,
					   Indonesia=Indonesia, Korea=Korea,
					   Russia=Russia, Australia=Australia, Germany=Germany )

	for (i in 1:length(conts)) {
		oo 	<- unlist(conts[i])
		
		for (j in 1:length(oo)) {
			inds 	<- which(places==oo[j])
			if (length(inds) > 0) {
				countries[inds] <- attributes(conts[i])$names
			}	 
		}
	}
	
	inds			<- which(countries=="-")
	countries[inds] 	<- places[inds]

	return (countries)

}


isProvinceOfChina	<- function( place ) {

	provinces <- c(	"Beijing","Tianjin","Hebei","Shanxi",
				"Mongolia","Liaoning","Jilin","Heilongjiang",
				"Shanghai","Jiangsu","Zhejiang","Anhui",
				"Fujian","Jiangxi","Shandong","Henan",
				"Hubei","Hunan","Guangdong",
				"Guangxi","Hainan","Chongqing","Sichuan",
				"Guizhou","Yunnan","Tibet","Xizang",
				"Shaanxi","Gansu","Qinghai","Ningxia",
				"Xinjiang","HongKong","Xianggang",
				"Macau","Aomen","Taiwan")
	provinces 	<- tolower(provinces)

	place2 	<- tolower(place)
	place2	<- gsub("_", "", place2)
	place2	<- gsub(" ", "", place2)
	place2	<- gsub("-", "", place2)
	
	if (length(place2) > 1) {
		inds		<- match(place2, provinces)
		jj		<- which(is.finite(inds))
		inds		<- array(FALSE, length(place2))
		inds[jj]	<- TRUE
		return( inds )
	} else {
		inds		<- which(provinces==place2)
		return ( length(inds) == 1 )
	}
}


# 2 Sept 2013 - updated
# 12 Feb 2014 - updated
getProvinceOfChina <- function( places ) {
	pc 		<- isProvinceOfChina( places )

	Anhui		<- c("Anhui","Chuzhou")
	Guizhou 	<- c("Guiyang")
	Zhejiang	<- c("Hangzhou","Huadong","Wenzhou")
	HongKong	<- c("HK","Hong_Kong")
	Jiangsu 	<- c(	"Jiawang","Nanjing","Qianzhou",
				"Taixing","Taizhou","Tongshan","Wuxi","Xiangshui","Xuzhou","Hongze",
				"Yangzhou","Shuanggou","Xigou","Xuyi","Suzhou")		# not sure about Shuanggou or Xigou
	Jiangxi 	<- c("Nanchang")
	Henan	  	<- c("Kaifeng","Qixian","Yongcheng")		# not sure about Qixian
	Yunnan  	<- c("Kuming","Kunming")
	Guangdong 	<- c("Shantou","Zhuhai","Shenzhen")
	Hebei		<- c("Shijiazhuang")
	Shanghai	<- c("Xianghai")
	Fujian	<- c("Zhangzhou")
	Shandong	<- c("Zibo","Rizhao")
	Shaanxi	<- c("Dawang")
	Guangxi	<- c("SanJiang","Sanjiang")				# not sure about SanJiang
	Heilongjiang<- c("ZhaLong")
	Hunan		<- c("Donting","Dongting")

	provinces	<- list( Anhui=Anhui, Guizhou=Guizhou, Zhejiang=Zhejiang, HongKong=HongKong,
				   Jiangsu=Jiangsu, Jiangxi=Jiangxi, Henan=Henan, Yunnan=Yunnan,
				   Guangdong=Guangdong, Hebei=Hebei, Shanghai=Shanghai, Fujian=Fujian,
				   Shandong=Shandong, Shaanxi=Shaanxi, 
				   Guangxi=Guangxi, Heilongjiang=Heilongjiang, Hunan=Hunan )

	provs			 <- array("-", length(places))
	provs[ which(pc) ] <- places[ which(pc) ]
	for (i in 1:length(provinces)) {
		oo 	<- unlist(provinces[i])
		
		for (j in 1:length(oo)) {
			inds 	<- which(places==oo[j])
			if (length(inds) > 0) {
				provs[inds] <- attributes(provinces[i])$names
			}	 
		}
	}
	
	return( provs )

}



