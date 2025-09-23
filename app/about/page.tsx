import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Dna, Users, Target, Heart, BookOpen, Lightbulb, Globe } from "lucide-react"

export default function AboutPage() {
  const team = [
    {
      name: "张三",
      role: "组织负责人",
      expertise: ["单细胞分析", "Python", "R"],
      description: "专注于单细胞转录组数据分析，有丰富的生物信息学项目经验。",
    },
    {
      name: "李四",
      role: "技术负责人",
      expertise: ["机器学习", "深度学习", "生物统计"],
      description: "致力于将机器学习方法应用于生物医学数据分析。",
    },
    {
      name: "王五",
      role: "内容负责人",
      expertise: ["基因组学", "转录组学", "数据可视化"],
      description: "专业从事基因组数据分析和生物信息学教学内容开发。",
    },
  ]

  const values = [
    {
      icon: BookOpen,
      title: "知识共享",
      description: "我们相信知识应该被自由分享，通过开放的学习资源帮助更多人掌握生物信息学技能。",
    },
    {
      icon: Users,
      title: "协作学习",
      description: "通过团队协作和同伴学习，我们能够更好地理解复杂的生物信息学概念和方法。",
    },
    {
      icon: Lightbulb,
      title: "创新实践",
      description: "鼓励创新思维，将最新的技术和方法应用到生物信息学研究和教学中。",
    },
    {
      icon: Globe,
      title: "开放包容",
      description: "欢迎来自不同背景的学生和研究者，共同构建多元化的学习社区。",
    },
  ]

  return (
    <div className="min-h-screen bg-background">
      <Navigation />

      {/* Hero Section */}
      <section className="py-20">
        <div className="container mx-auto px-4">
          <div className="max-w-4xl mx-auto text-center">
            <div className="flex justify-center mb-8">
              <div className="flex h-16 w-16 items-center justify-center rounded-2xl bg-primary text-primary-foreground">
                <Dna className="h-8 w-8" />
              </div>
            </div>
            <h1 className="text-4xl lg:text-5xl font-bold text-balance mb-6">关于我们</h1>
            <p className="text-xl text-muted-foreground text-pretty max-w-3xl mx-auto">
              上海生物信息中心学生组织成立于2020年，致力于为生物信息学学习者提供优质的学习资源、
              技术分享平台和学术交流机会。我们相信通过知识共享和协作学习， 能够推动生物信息学领域的发展和人才培养。
            </p>
          </div>
        </div>
      </section>

      {/* Mission & Vision */}
      <section className="py-16 bg-muted/50">
        <div className="container mx-auto px-4">
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            <Card>
              <CardHeader>
                <div className="flex items-center gap-3 mb-2">
                  <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-primary/10">
                    <Target className="h-5 w-5 text-primary" />
                  </div>
                  <CardTitle>我们的使命</CardTitle>
                </div>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground text-pretty">
                  通过提供高质量的学习资源、组织技术分享活动和建立学术交流平台，
                  帮助学生和研究者掌握生物信息学核心技能，培养创新思维， 推动生物信息学在生命科学研究中的应用和发展。
                </p>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <div className="flex items-center gap-3 mb-2">
                  <div className="flex h-10 w-10 items-center justify-center rounded-lg bg-primary/10">
                    <Heart className="h-5 w-5 text-primary" />
                  </div>
                  <CardTitle>我们的愿景</CardTitle>
                </div>
              </CardHeader>
              <CardContent>
                <p className="text-muted-foreground text-pretty">
                  成为国内领先的生物信息学学生社区，构建开放、包容、创新的学习环境，
                  培养具有国际视野和创新能力的生物信息学人才， 为推动生命科学研究的数字化转型贡献力量。
                </p>
              </CardContent>
            </Card>
          </div>
        </div>
      </section>

      {/* Core Values */}
      <section className="py-16">
        <div className="container mx-auto px-4">
          <div className="text-center mb-12">
            <h2 className="text-3xl font-bold mb-4">核心价值观</h2>
            <p className="text-muted-foreground max-w-2xl mx-auto">这些价值观指导着我们的行动，塑造着我们的社区文化</p>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
            {values.map((value, index) => {
              const Icon = value.icon
              return (
                <Card key={index}>
                  <CardContent className="p-6 text-center">
                    <div className="flex justify-center mb-4">
                      <div className="flex h-12 w-12 items-center justify-center rounded-lg bg-primary/10">
                        <Icon className="h-6 w-6 text-primary" />
                      </div>
                    </div>
                    <h3 className="text-lg font-semibold mb-3">{value.title}</h3>
                    <p className="text-sm text-muted-foreground text-pretty">{value.description}</p>
                  </CardContent>
                </Card>
              )
            })}
          </div>
        </div>
      </section>

      {/* Team */}
      <section className="py-16 bg-muted/50">
        <div className="container mx-auto px-4">
          <div className="text-center mb-12">
            <h2 className="text-3xl font-bold mb-4">核心团队</h2>
            <p className="text-muted-foreground max-w-2xl mx-auto">
              我们的团队由来自不同背景的生物信息学专家和学生组成
            </p>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            {team.map((member, index) => (
              <Card key={index}>
                <CardContent className="p-6">
                  <div className="flex items-center gap-4 mb-4">
                    <div className="h-12 w-12 rounded-full bg-primary/10 flex items-center justify-center">
                      <Users className="h-6 w-6 text-primary" />
                    </div>
                    <div>
                      <h3 className="text-lg font-semibold">{member.name}</h3>
                      <p className="text-sm text-muted-foreground">{member.role}</p>
                    </div>
                  </div>
                  <p className="text-sm text-muted-foreground mb-4 text-pretty">{member.description}</p>
                  <div className="flex flex-wrap gap-2">
                    {member.expertise.map((skill) => (
                      <Badge key={skill} variant="secondary" className="text-xs">
                        {skill}
                      </Badge>
                    ))}
                  </div>
                </CardContent>
              </Card>
            ))}
          </div>
        </div>
      </section>

      {/* Contact CTA */}
      <section className="py-20">
        <div className="container mx-auto px-4">
          <Card className="bg-primary text-primary-foreground">
            <CardContent className="p-12 text-center">
              <h2 className="text-3xl font-bold mb-4">加入我们</h2>
              <p className="text-lg opacity-90 mb-8 max-w-2xl mx-auto">
                如果你对生物信息学充满热情，希望与志同道合的伙伴一起学习和成长， 欢迎加入我们的团队。
              </p>
              <div className="flex flex-col sm:flex-row gap-4 justify-center">
                <a
                  href="mailto:contact@sbc.org.cn"
                  className="inline-flex items-center justify-center rounded-md bg-primary-foreground text-primary px-6 py-3 text-sm font-medium transition-colors hover:bg-primary-foreground/90"
                >
                  联系我们
                </a>
              </div>
            </CardContent>
          </Card>
        </div>
      </section>

      <Footer />
    </div>
  )
}
