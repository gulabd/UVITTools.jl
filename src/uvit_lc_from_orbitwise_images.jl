using FITSIO, Measurements, DataFrames, CSV


function uvit_lc_from_orbitwise_images(imagelist::String,srcds9regfile, bgdds9regfile, lcfile::String; mst_or_bjd::String="mst")
    imlist=open(imagelist)
    imfiles=readlines(imlist)
    tstart_mst=Array{Float64}(undef,length(imfiles))
    tstop_mst=Array{Float64}(undef, length(imfiles))
    time_meanmst=Array{Float64}(undef,length(imfiles))
    rat_with_error=Array{Float64}(undef,length(imfiles))
    cnt_rate=Array{Float64}(undef,length(imfiles))
    err_cnt_rate = Array{Float64}(undef,length(imfiles))
    flambda=Array{Float64}(undef,length(imfiles))
    err_flambda=Array{Float64}(undef,length(imfiles))
    err_mag=Array{Float64}(undef,length(imfiles))
    mag=Array{Float64}(undef,length(imfiles))

    for i in 1:length(imfiles)

        (time_meanmst[i], cnt_rate[i], err_cnt_rate[i]) = uvit_aphot(imfiles[i],srcds9regfile,bgdds9regfile,mst_or_bjd=mst_or_bjd)
        println(imfiles[i])
    end
  #  cnt_rate=Measurements.value(rat_with_error)
  #  err_cnt_rate=Measurements.uncertanty(rate_with_error)
  #  mag=Measurements.value(mag_with_error)
  #  err_mag=Measurements.value(mag_with_error)
  #  flambda=Measurements.value(flambda_with_error)
  #  err_flambda_mag=Measurements.value(flambda_with_error)
  #  time_sec=(time_meanbjd-time_meanbjd[1])*86400
    lc = DataFrame(mission_time=time_meanmst, cnt_rate=cnt_rate, err_cnt_rate=err_cnt_rate)
    # fluxlc=DataFrame(mission_time=time_meanmst, flambda=flambda, err_flambda=err_flambda)
    # maglc= DataFrame(mission_time=time_meanmst, mag=mag, err_mag=err_mag)
    #plot(time_meanbjd-time_meanbjd[1],rat,errorBars=ErrorBars(y=errs)
   # p=Axis(Plots.Linear(time_sec,cnt_rate,errorBars=ErrorBars(y=err_cnt_rate)),xlabel="Time (days)",ylabel=L"$\rm Counts~s^{-1}$")
   # save("lightcurve.svg",p)
    CSV.write("rate_" * lcfile,lc, delim=' ')
    # CSV.write("flux_" * lcfile,fluxlc, delim=' ')
  #  CSV.write("mag_lc" * lcfile,maglc, delim=' ')
    return time_meanmst,cnt_rate,err_cnt_rate
end
