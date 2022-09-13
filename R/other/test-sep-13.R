fun1 = function(a, b){
  print(eval(b))
}

b = expression(list(a, 1, 2))
fun1(a = "hello", b = b)

dep.var = "value"
data %>% filter(id != dependent.id) %>% group_by(time) %>% 
  summarise(average = mean(!!sym(dep.var), na.rm = TRUE)) %>% `$`(average)
data %>% filter(id != dependent.id) %>% group_by(time) %>% 
  summarise(average = mean(value, na.rm = TRUE)) %>% `$`(average)
