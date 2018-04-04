# Rate_control_schemes

The IEEE® 802.11™ standard supports dynamic rate control by adjusting the MCS value of each transmitted packet based on the underlying radio propagation channel. Maximizing link throughput, in a propagation channel that is time varying due to multipath fading or movement of the surrounding objects, requires dynamic variation of MCS. The IEEE 802.11 standard does not define any standardized rate control algorithm (RCA) for dynamically varying the modulation rate. The implementation of RCA is left open to the WLAN device manufacturers. This example uses a closed-loop rate control scheme. A recommended MCS for transmitting a packet is calculated at the receiver and is available at the transmitter without any feedback latency. In a real system this information would be conveyed through a control frame exchange. The MCS is adjusted for each subsequent packet in response to an evolving channel condition with noise power varying over time.


# Proposed scheme

We implemented link adaptation algorithm by calculating RSS_avg (average recived signal strength) and number of successfully transmitted packets to decide the MCS (Modulation and Coding scheme) to be used for next transmission. Algorith along with flow chart is explained in the report.



# references
1)https://nl.mathworks.com/help/wlan/examples/802-11-dynamic-rate-control-simulation.html
2) Q. Xia and M. Hamdi, "Smart sender: a practical rate adaptation algorithm for multirate IEEE 802.11 WLANs," in IEEE Transactions on Wireless Communications, vol. 7, no. 5, pp. 1764-1775, May 2008.
doi: 10.1109/TWC.2008.061047
3) J. P. Pavon and Sunghyun Choi, ”Link adaptation strategy for IEEE 802.11 WLAN via received signal strengthmeasurement,” Communications, 2003. ICC ’03. IEEE International Conference on, Anchorage, AK, 2003, pp. 1108-1113 vol.2. doi: 10.1109/ICC.2003.1204534
