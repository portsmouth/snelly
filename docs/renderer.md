
The renderer itself is a uni-directional pathtracer (with adjunct modes for ambient occlusion and normals rendering).

## Pathtracing
Relevent parameters are
```
  renderer.renderMode      (default 'pt')
  renderer.maxBounces      (default 4)
  renderer.maxMarchSteps 
  renderer.radianceClamp   (log scale)
  renderer.skyPower 
  renderer.skyTemperature 
  renderer.exposure 
  renderer.gamma 
  renderer.whitepoint 
  renderer.goalFPS 
```

### Lighting
For simplicity, the only light in the scene is a (non-HDRI) environment map. This can be specified via a URL to 
to a lat-long map, or otherwise will be taken to be a constant intensity sky. (See the envMap call).
In both cases, the sky spectrum is modulated by a blackbody emission spectrum with adjustable temperature.
