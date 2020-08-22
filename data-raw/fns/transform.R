transform_raw_aqol_tb_to_aqol6d_tb <- function(raw_aqol_tb){
  aqol6d_tb <- raw_aqol_tb %>%
    dplyr::mutate(d_agegroup=cut(d_age,breaks=c(11,17,30), labels=c("Age 12-17","Age 18-26"))) %>% # change
    dplyr::filter(!is.na(aqol6d_total_w)) %>%
    dplyr::mutate(round=factor(round,labels=c("Baseline","Follow-up"))) %>%
    dplyr::select(fkClientID, c_p_diag_s , s_centre, c_clinical_staging_s,
                  d_age , d_agegroup , d_gender,
                  d_sex_birth_s , d_sexual_ori_s ,
                  d_country_bir_s, d_ATSI,d_english_home, d_english_native,
                  d_relation_s,  d_studying_working, k6_total,
                  phq9_total, bads_total, gad7_total,oasis_total,
                  scared_total,
                  dplyr::contains("aqol6d"),c_sofas, round) %>%
    dplyr::mutate(Gender= factor(ifelse(d_gender =="Genderqueer/gender nonconforming/agender" |
                                          d_gender=="Transgender" , "Other", as.character(d_gender))))   %>%
    dplyr::mutate(Region=as.factor(ifelse(s_centre=="Canberra" |  s_centre=="Southport" | s_centre=="Knox", "Metro","Regional"))) %>%
    dplyr::mutate(CALD=factor(ifelse(d_country_bir_s=="Other" | d_english_home=="No" | d_english_native=="No" | d_ATSI=="Yes" , "Yes", "No" ))) %>%
    dplyr::rename(PHQ9=phq9_total,BADS=bads_total,GAD7=gad7_total,OASIS=oasis_total, SCARED=scared_total,K6=k6_total,SOFAS=c_sofas)
  Hmisc::label(aqol6d_tb$CALD)="Culturally and linguistically diverse (CALD) background"
  Hmisc::label(aqol6d_tb$d_agegroup)="Age group"
  return(aqol6d_tb)
}
