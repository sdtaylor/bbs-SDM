
fractions_skill_score = function(y_obs, y_pred){
  numerator   = mean((y_pred - y_obs)^2)
  denominator = sum(y_pred^2 + y_obs^2) / length(y_obs)
  return(1 - (numerator/denominator))
}

normalized_rmse = function(y_obs, y_pred){
  numerator = sqrt(mean((y_pred - y_obs)^2))
  denominator = max(y_obs) - min(y_obs)
  return(1-(numerator/denominator))
}