BEGIN{print "QSS,QSS_TT,DP,MQ,ReadPosRankSum,STVSB,SomaticEVS,TDP,TFDP,TSDP,TSUBDP,TAU,TCU,TGU,TTU,TDP,TFDP,TSDP,TSUBDP,TAU,TCU,TGU,TTU"}
$7~/PASS/ {
#!/^#/{
  # get values from the INFO field
  split($8,info,";")
  # we want to get         QSS, QSS_NT, DP, MQ, ReadPosRankSum, SNVSB, SomaticEVS from the info field
  # in the info_out array:  1    2      3   4    5               6       7
  for(i=1;i<=length(info);i++) {
    if (info[i] ~/^QSS=/) { split(info[i],qss,"="); out_info["qss"]=qss[2]; }
    if (info[i] ~/^QSS_NT=/)  { split(info[i],qss_nt,"="); out_info["qss_nt"]=qss_nt[2]; }
    if (info[i] ~/^DP=/) { split(info[i],dp,"="); out_info["dp"]=dp[2]; }
    if (info[i] ~/^MQ=/) { split(info[i],mq,"="); out_info["mq"]=mq[2]; }
    if (info[i] ~/RankSum=/) { split(info[i],ranksum,"="); out_info["ranksum"]=ranksum[2]; }
    if (info[i] ~/^SNVSB=/) { split(info[i],snvsb,"="); out_info["snvsb"]=snvsb[2]; }
    if (info[i] ~/^SomaticEVS=/) { split(info[i],sevs,"="); out_info["sevs"]=sevs[2]; }
  }

  # get genotype info - assuming NORMAL TUMOR order
  # DP:FDP:SDP:SUBDP:AU:CU:GU:TU
  # FORMAT  $9
  # NORMAL $10
  # TUMOR  $11

  # simple sub() replacement is not enough, as we have list of values sometimes
  # for these list now we will get the first value only
  split($10,normal,":")
  # DP, FDP, SDP  and SUBDP are single
  ngt["DP"] = normal[1]
  ngt["FDP"] = normal[2]
  ngt["SDP"] = normal[3]
  ngt["SUBDP"] = normal[4]
  # get the first for these
  split(normal[1],au,","); ngt["AU"]=au[1];
  split(normal[1],cu,","); ngt["CU"]=cu[1];
  split(normal[1],gu,","); ngt["GU"]=gu[1];
  split(normal[1],tu,","); ngt["TU"]=tu[1];
  normal_gt = ngt["DP"]","ngt["FDP"]","ngt["SDP"]","ngt["SUBDP"]","ngt["AU"]","ngt["CU"]","ngt["GU"]","ngt["TU"]

  # ditto for tumor
  split($11,tumor,":")
  # DP, FDP, SDP  and SUBDP are single
  tgt["DP"] = tumor[1]
  tgt["FDP"] = tumor[2]
  tgt["SDP"] = tumor[3]
  tgt["SUBDP"] = tumor[4]
  # get the first for these
  split(tumor[1],au,","); tgt["AU"]=au[1];
  split(tumor[1],cu,","); tgt["CU"]=cu[1];
  split(tumor[1],gu,","); tgt["GU"]=gu[1];
  split(tumor[1],tu,","); tgt["TU"]=tu[1];
  tumor_gt = tgt["DP"]","tgt["FDP"]","tgt["SDP"]","tgt["SUBDP"]","tgt["AU"]","tgt["CU"]","tgt["GU"]","tgt["TU"]

  info_string = out_info["qss"]","out_info["qss_nt"]","out_info["dp"]","out_info["mq"]","out_info["ranksum"]","out_info["snvsb"]","out_info["sevs"]


  print info_string","normal_gt","tumor_gt 
}
