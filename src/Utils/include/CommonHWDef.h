#ifndef COMMONHWDEF_H
#define COMMONHWDEF_H

#include "Definitions.h"
#define constLogPrfHeader  "\xba\xdc\xfe\x21"
#define constDspCardIp  "192.168.100.3"
#define constDspCardPort  1025
#define constExtraFooterSamples  (sizeof(PrfPacketFooter) / sizeof(int16_t)) // bytes

#define CONST_MIN_PRF_READ 2*5000
#define CONST_RAW_BUFFER_SIZE 20*1024*1024
#define ConstStatusReaderTimeout  2000

#define CONST_CORE_SETTING_FILE_NAME "CoreSetting.db"

struct LongRangePrfPack
{
	int seqNumber;
	int prfLen;
	char prfData[100*1024];
};

struct MotionDataInput
{
// 	MotionParametes motionParams;
	double st_m, en_m, numPulses, Z, powr, powa, powz, fac;
};
struct MotionDataOutput
{
	double vel,a,r,z;
};

enum InputType
{
    InputType_Udp = 0,
	InputType_File ,
	InputType_TCP ,
	InputType_TCP_FILE ,
	InputType_TCPV3 ,
	InputType_TCP_FILEV3 ,
	
};

//#define LogRecorder_Last_Browse_Dir  "LastBrowseDir"
//#define Param_LocalPort			 "LocalPort"
#define Param_LocalPortUdp			 "LocalPortUdp"
#define Param_InputFile			 "InputFile"
#define Param_InputType			 "InputTypeNetOrFile"
#define Param_AutoReload			 "AutoReload"
#define Param_MovingPulses_SAR	 "MovingPulsesSar"
// #define Param_MovingPulses_SAR_rate	 "MovingPulsesSar_rate"
#define Param_MovingPulses_GMTI	 "MovingPulsesGmti"
#define Param_MovingPulses_WGMTI	 "wgmti_moving_pulse"
#define Param_MovingPulses_MMTI	 "mmti_moving_pulse"
#define Param_MovingPulses_ISAR	 "ISAR_moving_pulse"
#define Param_FirstPulses_SAR	 "FirstPulsesSar"
#define Param_FirstPulses_GMTI	 "FirstPulsesGmti"
#define Param_FirstPulses_WGMTI	 "wgmti_total_pulse"
#define Param_FirstPulses_MMTI	 "mmti_total_pulse"
#define Param_FirstPulses_ISAR	 "isar_total_pulse"
#define Param_NStart_SAR			 "NStartSar"
#define Param_NStop_SAR			 "NStopSar"
#define Param_NStart_GMTI		 "NStartGmti"
#define Param_NStop_GMTI			 "NStopGmti"
#define Param_GPS_TIME           "CheckGPSTime"
#define Param_NStart_WGMTI		 "NStartWGmti"
#define Param_NStop_WGMTI		 "NStopWGmti"
#define Param_SweepBw_WGMTI		 "SweepBwWGmti"
#define Param_DataBw_WGMTI		 "DataBwWGmti"
#define Param_FirstPulses_SAR_Manual	 "FirstPulsesSarManual"

#define Param_CardTcpIp			 "CardTcpIp"
#define Param_CardTcpPort		 "CardTcpPort"
#define Param_LogDir				 "LogDir"
//#define Param_SaveFlightPath		 "saveFlightPath"
#define Param_DelayBetweenReadingUs		 "delayReadingPulseUs"
#define Param_GmtiLimitPulseLen		 "gmtiLimitPulseLen"
//#define Param_GmtisLimitPulseLen		 "gmtisLimitPulseLen"
#define Param_ChecksumCheck		 "checksumCheck"
//#define Param_GmtiAutoMode		 "gmtiauto"
//#define Param_CircularAutoMode   "circularauto"
//#define Param_GmtisAutoMode		 "gmtisauto"
//#define Param_MmtisAutoMode		 "mmtisauto"
//static const char* Param_GmtiAutoModeRefLoc		 "gmtiautorefloc"
//static const char* Param_ImagingLocA		 "imagingloca"
//static const char* Param_ImagingLocB		 "imaginglocb"
//static const char* Param_ImagingLocC		 "imaginglocc"
//static const char* Param_ImagingH		 "imagingh"
//#define Param_ImagingAutoMode		 "imagingauto"
#define Param_OffsetAltitude		 "offsetaltitude"
//#define Param_DictGmtiScenario		 "gmtiscenariolist"
//#define Param_DictGmtisScenario		 "gmtisscenariolist"
#define Param_DICT_MMTIS_SCENARIO	 "mmtisscenariolist"
#define Param_GmtiScenarioName		 "gmtiscenarioname"
#define Param_mti_search_scenario_name		 "gmtisscenarioname"
// #define Param_MMTIS_SCENARIO_NAME		 "mmtisscenarioname"
//#define Param_DictImagingScenario		 "imagingscenariolist"
#define Param_ImagingScenarioName		 "imagingscenarioname"
#define Param_ISarImagingScenarioName		 "isarimagingscenarioname"

#define Param_AutoModeType          "automode_type"
#define Param_ReceiverBandwidth     "receiver_bandwidth"

#define Param_CircularScenarioName		 "circularscenarioname"
#define Param_DictCircularScenario		 "circularscenariolist"
#define Param_UseFixedYaw               "usefixedyaw"
#define Param_FixedYaw                  "fixedyaw"
#define Param_ServoDown                  "servodown"
#define PresumGpsFrq  "GpsFrq"
#define ImagingAckPack "imagingackpack2"
#define GmtiAckPack "gmtiackpack2"
#define GmtisAckPack "gmtisackpack2"
#define WGmtisAckPack "wgmtisackpack2"
#define MmtisAckPack "mmtisackpack2"
#define Param_ScenarioOverride "scenario_override"
#define Param_CoreLogsAddress "CoreLogeAddress"
#define Param_AutoRecordRaw "auto_record_raw"
#define Param_AutoScenario "auto_Scenario_enable"
#define Param_PORT_UDP_AIS			 "ais_port"
#define Param_GPS_Lock           "CheckGPSLock"
#define Param_PAN_CHECK           "CheckPan"


#define Presum_HeightMapFile2 "HeightMapFile"
#define ImageHeightMode "image_height_mode"
#define IRANHEIGHTADDRESS "/home/sarjahad/map/iranheight.tif"

//static const char* Param_FixedAreaAutoTilt		= "imagingauto";

#endif // COMMONHWDEF_H
