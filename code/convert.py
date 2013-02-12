import string
f=open('train-hoghof.csv','w')
f.write('hog1, hog2, hog3, hog4, hog5, hog6, hog7, hog8, hog9, hog10, hog11, hog12, hog13, hog14, hog15, hog16, hog17, hog18, hog19, hog20, hog21, hog22, hog23, hog24, hog25, hog26, hog27, hog28, hog29, hog30, hog31, hog32, hog33, hog34, hog35, hog36, hog37, hog38, hog39, hog40, hog41, hog42, hog43, hog44, hog45, hog46, hog47, hog48, hog49, hog50, hog51, hog52, hog53, hog54, hog55, hog56, hog57, hog58, hog59, hog60, hog61, hog62, hog63, hog64, hog65, hog66, hog67, hog68, hog69, hog70, hog71, hog72, hof1, hof2, hof3, hof4, hof5, hof6, hof7, hof8, hof9, hof10, hof11, hof12, hof13, hof14, hof15, hof16, hof17, hof18, hof19, hof20, hof21, hof22, hof23, hof24, hof25, hof26, hof27, hof28, hof29, hof30, hof31, hof32, hof33, hof34, hof35, hof36, hof37, hof38, hof39, hof40, hof41, hof42, hof43, hof44, hof45, hof46, hof47, hof48, hof49, hof50, hof51, hof52, hof53, hof54, hof55, hof56, hof57, hof58, hof59, hof60, hof61, hof62, hof63, hof64, hof65, hof66, hof67, hof68, hof69, hof70, hof71, hof72, hof73, hof74, hof75, hof76, hof77, hof78, hof79, hof80, hof81, hof82, hof83, hof84, hof85, hof86, hof87, hof88, hof89, hof90\n')
for line in open('train-hoghof.txt','r'):
	if line[0]!='#':
		for x in line.split()[9:-1]:
			f.write(str(int(round(float(x)*255))) + ',')
		f.write(str(int(round(float(x)*255))) + '\n')
