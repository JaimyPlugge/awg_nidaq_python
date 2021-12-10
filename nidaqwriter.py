"""
This file contains a Writer class that uses nidaqmx
to communicate with NI DAQs that are capable of
outputting a voltage. This module was especially tested 
on a NI myDAQ, but should also work on other output DAQs.

@author: JaimyPlugge
"""
import nidaqmx
import numpy as np
from nidaqmx import stream_writers

class Writer: 
    def __init__(self, chan_name1, chan_name2, sample_rate):
        self.sample_rate = sample_rate
        self.chan_name1 = chan_name1
        self.chan_name2 = chan_name2
        self.task = nidaqmx.Task()
        self.task.ao_channels.add_ao_voltage_chan(self.chan_name1)
        self.multichan = False
        if len(self.chan_name2) > 0:
            self.multichan = True
            self.task.ao_channels.add_ao_voltage_chan(self.chan_name2)


    def changetask(self, new_chan_name1, new_chan_name2):
        self.stopfunc()
        self.chan_name1 = new_chan_name1
        self.chan_name2 = new_chan_name2
        self.task = nidaqmx.Task()
        self.task.ao_channels.add_ao_voltage_chan(self.chan_name1)
        self.multichan = False
        if len(self.chan_name2) > 0:
            self.multichan = True
            self.task.ao_channels.add_ao_voltage_chan(self.chan_name2)

        
    def outputcontinuously(self, waveform):
        self.task.timing.cfg_samp_clk_timing(rate= self.sample_rate, sample_mode= nidaqmx.constants.AcquisitionType.CONTINUOUS, samps_per_chan= 10)
        if self.multichan:
            test_Writer = stream_writers.AnalogMultiChannelWriter(self.task.out_stream, auto_start=False)
            test_Writer.write_many_sample(waveform, timeout=nidaqmx.constants.WAIT_INFINITELY)
        else:    
            test_Writer = stream_writers.AnalogSingleChannelWriter(self.task.out_stream, auto_start=False)
            test_Writer.write_many_sample(waveform[0,:], timeout=nidaqmx.constants.WAIT_INFINITELY)
        self.task.start()

    def singleoutput(self, samples):
        self.task.timing.cfg_samp_clk_timing(rate= self.sample_rate, sample_mode= nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan= samples.shape[1])
        if self.multichan:
            test_Writer = stream_writers.AnalogMultiChannelWriter(self.task.out_stream, auto_start=False)
            test_Writer.write_many_sample(samples)
        else:    
            test_Writer = stream_writers.AnalogSingleChannelWriter(self.task.out_stream, auto_start=False)
            test_Writer.write_many_sample(samples[0,:])
        self.task.start()


    def pausefunc(self):
        self.task.stop()

    def stopfunc(self):
        self.pausefunc()
        self.outputcontinuously(np.zeros((2, 10), dtype=float))
        self.task.close()