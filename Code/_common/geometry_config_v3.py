#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
geometry_config_v3.py
=====================

秦皇岛 barrier–lagoon DAS 几何配置。

这份配置把两类信息严格区分：

A. “确定/强约束”的几何信息
--------------------------
来自通道号与 2 m 通道间距：
- Lagoon1: CH32–101
- Lagoon2: CH102–163
- Lagoon3: CH164–257
- Connector: CH257–285–310–330
- Coast: CH330–430
- Coast extension: CH430–450
- Channel spacing = 2 m

因此：
- CH32–101: 69 intervals = 138 m
- CH101–163: 62 intervals = 124 m
- CH163–257: 94 intervals = 188 m
- CH257–285: 28 intervals = 56 m
- CH285–310: 25 intervals = 50 m
- CH310–330: 20 intervals = 40 m
- CH330–430: 100 intervals = 200 m
- CH430–450: 20 intervals = 40 m

Lagoon 三角形三边 138/124/188 m 会在代码中严格满足，
相应三个内角由余弦定理确定，约 91.5°, 47.2°, 41.3°。

B. “地图近似”的整体方位信息
--------------------------
下面几个角度来自你提供的布设图，用于把严格长度几何旋转/摆放到
接近真实地图的方向。它们不是 RTK/GIS 真值，因此代码会：
1. 在输出中明确标记为 map-derived approximate azimuth；
2. 提供 11_几何不确定性敏感性.py，对 ±5° 方位不确定性做 Monte Carlo。

角度定义：
- 数学平面角，0° = east, 90° = north
- 方向只在 0–180° 中有意义，因为 DAS 光纤轴向无正反区别

如果未来得到 RTK/GIS 拐点坐标，只需改这份配置或直接替换几何表。
"""

CHANNEL_SPACING_M = 2.0

# analysis segments
SITE_CHANNEL_RANGES = {
    "Lagoon1": (32, 101),
    "Lagoon2": (102, 163),
    "Lagoon3": (164, 257),
    "Coast":   (330, 430),
}

# exact path breakpoints from deployment
PATH_BREAKS = [32, 101, 163, 257, 285, 310, 330, 430, 450]

# exact segment lengths from channel difference * 2 m
EXACT_SEGMENT_LENGTHS_M = {
    (32, 101): 138.0,
    (101, 163): 124.0,
    (163, 257): 188.0,
    (257, 285): 56.0,
    (285, 310): 50.0,
    (310, 330): 40.0,
    (330, 430): 200.0,
    (430, 450): 40.0,
}

# Lagoon triangle overall rotation:
# CH101 -> CH163 is close to east-west and rises slightly eastward on the map.
TOP_EDGE_AZIMUTH_DEG_FROM_EAST = 6.0

# Map-derived approximate directional bearings for connector/coast segments.
# These preserve exact lengths while approximating the layout figure.
MAP_SEGMENT_AZIMUTHS_DEG = {
    (257, 285): -25.0,   # southeast from lagoon lower vertex
    (285, 310): -136.0,  # southwest
    (310, 330): -48.0,   # southeast
    (330, 430): 63.0,    # northeast, Coast low-tide-line segment
    (430, 450): 158.0,   # west-northwest short top extension
}

# The map shows CH32 and CH257 returning to essentially the same lower lagoon vertex.
# With only the figure we cannot justify a precise few-metre separation, so nominal geometry
# treats them as coincident. Sensitivity analysis can perturb this if needed later.
LOWER_VERTEX_GAP_M = 0.0

# uncertainty used by optional sensitivity analysis
MAP_AZIMUTH_UNCERTAINTY_DEG = 5.0
GLOBAL_ROTATION_UNCERTAINTY_DEG = 3.0

# core representative channels
REPRESENTATIVE_CHANNELS = {
    "Coast": 413,
    "Lagoon1": 67,
    "Lagoon2": 130,
    "Lagoon3": 186,
}
