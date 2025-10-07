import { Navigation } from "@/components/navigation"
import { Footer } from "@/components/footer"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Button } from "@/components/ui/button"
import { Dna, Users, Target, Heart, BookOpen, Lightbulb, Globe, Mail, MessageCircle, Calendar } from "lucide-react"
import Image from "next/image"
import { getPosts, getAnnouncements, getLearnResources, getMeetings } from "@/lib/server-markdown-loader"

export default async function AboutPage() {
  const posts = getPosts()
  const announcements = getAnnouncements()
  const resources = getLearnResources()
  const meetings = getMeetings()

  const stats = {
    posts: posts.length,
    announcements: announcements.length,
    resources: resources.length,
    meetings: meetings.length
  }
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
      title: "开放共享",
      description: "知识应该被自由分享，通过开放的学习资源和成果展示帮助更多人了解生物信息学。",
      emoji: "📚"
    },
    {
      icon: Users,
      title: "合作共赢",
      description: "通过团队协作和同伴交流，我们能够更好地解决复杂的生物信息学问题。",
      emoji: "🤝"
    },
    {
      icon: Lightbulb,
      title: "创新思维",
      description: "鼓励创新思考，勇于挑战传统，追求突破性进展。",
      emoji: "💡"
    },
    {
      icon: Globe,
      title: "学术自由",
      description: "欢迎来自不同背景的学生和研究者，构建多元化的学术社区。",
      emoji: "🌍"
    },
  ]

  const activities = [
    {
      icon: Calendar,
      title: "学术研讨会",
      description: "每周举行不同主题的学术讨论，分享最新研究成果",
      emoji: "🎯"
    },
    {
      icon: MessageCircle,
      title: "技术分享",
      description: "成员定期分享编程技巧和数据分析经验",
      emoji: "💻"
    },
    {
      icon: Users,
      title: "合作项目",
      description: "组队参与科研项目，在实践中提升能力",
      emoji: "🚀"
    },
  ]

  return (
    <div className="min-h-screen bg-background">
      <Navigation />

      {/* Hero Section - Tech Style */}
      <section className="relative min-h-screen flex items-center overflow-hidden tech-gradient">
        {/* Animated Background Grid */}
        <div className="absolute inset-0 tech-grid opacity-30"></div>

        {/* Floating Elements */}
        <div className="absolute top-20 left-10 w-20 h-20 bg-primary/10 rounded-full blur-xl animate-blob"></div>
        <div className="absolute bottom-20 right-10 w-32 h-32 bg-accent/10 rounded-full blur-xl animate-blob animation-delay-2000"></div>
        <div className="absolute top-1/2 left-1/4 w-16 h-16 bg-primary/20 rounded-full blur-lg animate-blob animation-delay-4000"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="max-w-7xl mx-auto">
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-16 items-center">
              {/* Left Content */}
              <div className="space-y-8">
                {/* Logo and Title */}
                <div className="space-y-6">
                  <div className="flex items-center gap-4">
                    <div className="relative h-24 w-24 pulse-glow rounded-2xl bg-primary/10 p-2" style={{ position: 'relative' }}>
                      <Image
                        src="/logo-nav.png"
                        alt="SBC Logo"
                        fill
                        className="object-contain"
                      />
                    </div>
                    <div>
                      <div className="text-sm font-mono text-primary/60 mb-2">SINCE 2025</div>
                      <h1 className="text-5xl lg:text-7xl font-bold text-balance">
                        <span className="block">ABOUT</span>
                        <span className="block text-primary neon-glow">SBC</span>
                      </h1>
                    </div>
                  </div>

                  <div className="space-y-4">
                    <h2 className="text-2xl lg:text-3xl font-bold text-foreground/90">
                      上海生物信息学中心
                      <span className="block text-primary">学生学术联盟</span>
                    </h2>

                    <p className="text-lg text-foreground/70 leading-relaxed">
                      我们是一个由热爱生物信息学的学生自发组织的
                      <span className="text-primary font-semibold">前沿科技学术社区</span>，
                      致力于促进学术交流、推动技术创新和培养未来生物信息学人才。
                    </p>
                  </div>

                  {/* Tech Tags */}
                  <div className="flex flex-wrap gap-3">
                    {[
                      "🧬 生物信息学",
                      "👥 学术交流",
                      "💡 技术创新",
                      "🎓 人才培养",
                      "🔬 科研合作",
                      "🚀 前沿探索"
                    ].map((tag, index) => (
                      <span
                        key={index}
                        className="px-4 py-2 rounded-full bg-primary/10 text-primary border border-primary/20
                               hover:bg-primary/20 transition-all duration-300 hover:scale-105
                               font-mono text-sm cyber-border"
                      >
                        {tag}
                      </span>
                    ))}
                  </div>

                  {/* Stats */}
                  <div className="grid grid-cols-3 gap-6 pt-6">
                    {[
                      { number: `${stats.activeMembers}+`, label: "活跃成员" },
                      { number: `${stats.techShares}+`, label: "技术分享" },
                      { number: `${stats.collaborationProjects}+`, label: "合作项目" }
                    ].map((stat, index) => (
                      <div key={index} className="text-center">
                        <div className="text-3xl font-bold text-primary mb-1">{stat.number}</div>
                        <div className="text-sm text-foreground/60 font-mono">{stat.label}</div>
                      </div>
                    ))}
                  </div>
                </div>
              </div>

              {/* Right Content - Animated Characters */}
              <div className="relative">
                {/* Main Characters Container */}
                <div className="relative flex justify-center items-center gap-8">
                  {/* Character 1 */}
                  <div className="relative w-64 h-64 animate-float" style={{ position: 'relative' }}>
                    <div className="absolute inset-0 bg-gradient-to-br from-primary/20 to-accent/20 rounded-full blur-2xl"></div>
                    <Image
                      src="/images/上海生物信息中心卡通形象.png"
                      alt="SBC Tech Mascot 1"
                      fill
                      className="object-contain relative z-10 drop-shadow-2xl filter hover:scale-110 transition-transform duration-300"
                    />
                  </div>

                  {/* Character 2 */}
                  <div className="relative w-64 h-64 animate-float-delayed" style={{ position: 'relative' }}>
                    <div className="absolute inset-0 bg-gradient-to-br from-accent/20 to-primary/20 rounded-full blur-2xl"></div>
                    <Image
                      src="/images/上海生物信息中心卡通形象v2.png"
                      alt="SBC Tech Mascot 2"
                      fill
                      className="object-contain relative z-10 drop-shadow-2xl filter hover:scale-110 transition-transform duration-300"
                    />
                  </div>
                </div>

                {/* Floating Tech Elements */}
                <div className="absolute -top-10 -right-10 w-20 h-20 bg-primary/10 rounded-lg rotate-45 animate-pulse"></div>
                <div className="absolute -bottom-10 -left-10 w-16 h-16 bg-accent/10 rounded-lg rotate-12 animate-pulse animation-delay-2000"></div>
              </div>
            </div>
          </div>
        </div>

        {/* Scroll Indicator */}
        <div className="absolute bottom-8 left-1/2 transform -translate-x-1/2 text-center">
          <div className="w-6 h-10 border-2 border-primary/30 rounded-full mx-auto mb-2">
            <div className="w-1 h-3 bg-primary/60 rounded-full mx-auto mt-2 animate-bounce"></div>
          </div>
          <p className="text-xs text-primary/60 font-mono">SCROLL TO EXPLORE</p>
        </div>
      </section>

      {/* Mission & Vision - Tech Style */}
      <section className="py-24 relative bg-gradient-to-br from-background via-primary/5 to-background">
        {/* Tech Background */}
        <div className="absolute inset-0 tech-grid opacity-20"></div>
        <div className="container mx-auto px-4 relative z-10">
          <div className="max-w-6xl mx-auto">
            <div className="text-center mb-16">
              <h2 className="text-4xl font-bold mb-4">
                <span className="text-foreground/90">OUR</span>
                <span className="text-primary neon-glow"> MISSION</span>
              </h2>
              <p className="text-lg text-foreground/60 font-mono">驱动创新 · 连接未来 · 培养人才</p>
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-12">
              {/* Mission */}
              <div className="group">
                <div className="glass-effect rounded-2xl p-8 border border-primary/20
                            hover:border-primary/40 transition-all duration-500
                            hover:shadow-2xl hover:shadow-primary/10">
                  <div className="flex items-center gap-4 mb-6">
                    <div className="flex h-16 w-16 items-center justify-center rounded-xl bg-primary/10
                                group-hover:bg-primary/20 transition-all duration-300">
                      <Target className="h-8 w-8 text-primary" />
                    </div>
                    <h3 className="text-2xl font-bold text-primary">MISSION</h3>
                  </div>
                  <div className="space-y-4">
                    <p className="text-foreground/80 leading-relaxed text-lg">
                      搭建学生间学术交流的<span className="text-primary font-semibold">桥梁</span>，
                      促进知识分享和思想碰撞
                    </p>
                    <p className="text-foreground/70 leading-relaxed">
                      培养具备创新能力的生物信息学人才，
                      推动生物信息学技术在生命科学领域的应用和发展
                    </p>
                  </div>
                </div>
              </div>

              {/* Vision */}
              <div className="group">
                <div className="glass-effect rounded-2xl p-8 border border-accent/20
                            hover:border-accent/40 transition-all duration-500
                            hover:shadow-2xl hover:shadow-accent/10">
                  <div className="flex items-center gap-4 mb-6">
                    <div className="flex h-16 w-16 items-center justify-center rounded-xl bg-accent/10
                                group-hover:bg-accent/20 transition-all duration-300">
                      <Heart className="h-8 w-8 text-accent" />
                    </div>
                    <h3 className="text-2xl font-bold text-accent">VISION</h3>
                  </div>
                  <div className="space-y-4">
                    <p className="text-foreground/80 leading-relaxed text-lg">
                      成为具有<span className="text-accent font-semibold">影响力</span>的学生学术联盟，
                      构建开放包容的学术社区
                    </p>
                    <p className="text-foreground/70 leading-relaxed">
                      培养具有创新精神和实践能力的生物信息学人才，
                      为推动生命科学研究和技术创新贡献力量
                    </p>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </section>

      {/* Core Values - Tech Style */}
      <section className="py-24 relative">
        {/* Animated Background */}
        <div className="absolute inset-0 tech-gradient opacity-30"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="max-w-6xl mx-auto">
            <div className="text-center mb-20">
              <h2 className="text-5xl font-bold mb-6">
                <span className="text-foreground/90">CORE</span>
                <span className="text-primary neon-glow"> VALUES</span>
              </h2>
              <p className="text-lg text-foreground/60 font-mono max-w-2xl mx-auto">
                这些价值观指导着我们的行动，塑造着我们独特的科技社区文化
              </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-8">
              {values.map((value, index) => {
                const Icon = value.icon
                return (
                  <div key={index} className="group relative">
                    {/* Glass Card */}
                    <div className="glass-effect rounded-2xl p-8 border border-primary/10
                                hover:border-primary/30 transition-all duration-500
                                hover:shadow-2xl hover:shadow-primary/10
                                hover:-translate-y-2 transform">
                      {/* Icon and Emoji */}
                      <div className="text-center space-y-6">
                        <div className="relative">
                          <div className="text-6xl mb-4 group-hover:scale-110 transition-transform duration-300">
                            {value.emoji}
                          </div>
                          <div className="absolute -top-2 -right-2 w-8 h-8 bg-primary/20 rounded-full blur-sm group-hover:bg-primary/40 transition-colors"></div>
                        </div>

                        <div className="flex justify-center">
                          <div className="flex h-12 w-12 items-center justify-center rounded-xl bg-primary/10
                                      group-hover:bg-primary/20 transition-all duration-300">
                            <Icon className="h-6 w-6 text-primary" />
                          </div>
                        </div>

                        <div className="space-y-3">
                          <h3 className="text-xl font-bold text-primary group-hover:text-primary/80 transition-colors">
                            {value.title}
                          </h3>
                          <p className="text-sm text-foreground/70 leading-relaxed">
                            {value.description}
                          </p>
                        </div>
                      </div>
                    </div>

                    {/* Hover Effect Border */}
                    <div className="absolute inset-0 rounded-2xl border-2 border-transparent
                                group-hover:border-primary/30 transition-all duration-500
                                pointer-events-none"></div>
                  </div>
                )
              })}
            </div>
          </div>
        </div>
      </section>

      {/* Activities - Tech Style */}
      <section className="py-24 relative bg-gradient-to-br from-background via-accent/5 to-background">
        {/* Animated Background */}
        <div className="absolute inset-0 tech-grid opacity-20"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="max-w-6xl mx-auto">
            <div className="text-center mb-20">
              <h2 className="text-5xl font-bold mb-6">
                <span className="text-foreground/90">OUR</span>
                <span className="text-accent neon-glow"> ACTIVITIES</span>
              </h2>
              <p className="text-lg text-foreground/60 font-mono max-w-2xl mx-auto">
                定期举办各种前沿学术和科技活动，丰富成员的学习体验
              </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
              {activities.map((activity, index) => {
                const Icon = activity.icon
                return (
                  <div key={index} className="group">
                    <div className="glass-effect rounded-2xl p-8 border border-accent/20
                                hover:border-accent/40 transition-all duration-500
                                hover:shadow-2xl hover:shadow-accent/10
                                hover:-translate-y-2 transform h-full">
                      <div className="text-center space-y-6">
                        {/* Emoji with glow effect */}
                        <div className="relative">
                          <div className="text-5xl group-hover:scale-110 transition-transform duration-300">
                            {activity.emoji}
                          </div>
                          <div className="absolute inset-0 bg-accent/20 rounded-full blur-xl opacity-0 group-hover:opacity-100 transition-opacity"></div>
                        </div>

                        {/* Icon */}
                        <div className="flex justify-center">
                          <div className="flex h-14 w-14 items-center justify-center rounded-xl bg-accent/10
                                      group-hover:bg-accent/20 transition-all duration-300">
                            <Icon className="h-7 w-7 text-accent" />
                          </div>
                        </div>

                        {/* Content */}
                        <div className="space-y-4">
                          <h3 className="text-xl font-bold text-accent group-hover:text-accent/80 transition-colors">
                            {activity.title}
                          </h3>
                          <p className="text-sm text-foreground/70 leading-relaxed">
                            {activity.description}
                          </p>
                        </div>

                        {/* Tech Decorative Element */}
                        <div className="flex justify-center">
                          <div className="w-16 h-1 bg-gradient-to-r from-accent/20 to-accent/60 rounded-full"></div>
                        </div>
                      </div>
                    </div>
                  </div>
                )
              })}
            </div>
          </div>
        </div>
      </section>

      {/* Team - Tech Style */}
      <section className="py-24 relative">
        {/* Animated Background */}
        <div className="absolute inset-0 tech-gradient opacity-20"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="max-w-6xl mx-auto">
            <div className="text-center mb-20">
              <h2 className="text-5xl font-bold mb-6">
                <span className="text-foreground/90">CORE</span>
                <span className="text-primary neon-glow"> TEAM</span>
              </h2>
              <p className="text-lg text-foreground/60 font-mono max-w-2xl mx-auto">
                联盟由热爱生物信息学的学生自发组织和维护
              </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
              {team.map((member, index) => (
                <div key={index} className="group">
                  <div className="glass-effect rounded-2xl p-8 border border-primary/10
                              hover:border-primary/30 transition-all duration-500
                              hover:shadow-2xl hover:shadow-primary/10
                              hover:-translate-y-2 transform h-full">
                    <div className="text-center space-y-6">
                      {/* Avatar with glow effect */}
                      <div className="relative mx-auto">
                        <div className="relative w-24 h-24 mx-auto" style={{ position: 'relative' }}>
                          <div className="absolute inset-0 bg-primary/20 rounded-full blur-xl opacity-0 group-hover:opacity-100 transition-opacity"></div>
                          <Image
                            src={`/images/上海生物信息中心表情包${index + 1 === 1 ? 'v2' : index + 1 === 2 ? '' : 'v2'}.png`}
                            alt={member.name}
                            fill
                            className="object-contain relative z-10 rounded-full group-hover:scale-110 transition-transform duration-300"
                          />
                        </div>
                        {/* Tech ring around avatar */}
                        <div className="absolute inset-0 w-32 h-32 mx-auto border-2 border-primary/20 rounded-full group-hover:border-primary/40 transition-colors"></div>
                      </div>

                      {/* Member Info */}
                      <div className="space-y-4">
                        <div>
                          <h3 className="text-xl font-bold text-primary group-hover:text-primary/80 transition-colors">
                            {member.name}
                          </h3>
                          <p className="text-sm text-foreground/60 font-mono uppercase tracking-wide">
                            {member.role}
                          </p>
                        </div>

                        <p className="text-sm text-foreground/70 leading-relaxed">
                          {member.description}
                        </p>

                        {/* Skills */}
                        <div className="flex flex-wrap gap-2 justify-center">
                          {member.expertise.map((skill) => (
                            <span
                              key={skill}
                              className="px-3 py-1 rounded-full bg-primary/10 text-primary text-xs font-mono
                                     hover:bg-primary/20 transition-colors border border-primary/20"
                            >
                              {skill}
                            </span>
                          ))}
                        </div>
                      </div>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        </div>
      </section>

      {/* Contact CTA - Tech Style */}
      <section className="relative py-32 overflow-hidden">
        {/* Animated Background */}
        <div className="absolute inset-0 tech-gradient opacity-40"></div>
        <div className="absolute inset-0 tech-grid opacity-30"></div>

        {/* Floating Elements */}
        <div className="absolute top-10 left-10 w-24 h-24 bg-primary/20 rounded-full blur-2xl animate-blob"></div>
        <div className="absolute bottom-10 right-10 w-32 h-32 bg-accent/20 rounded-full blur-2xl animate-blob animation-delay-2000"></div>
        <div className="absolute top-1/2 left-1/3 w-20 h-20 bg-primary/30 rounded-full blur-xl animate-blob animation-delay-4000"></div>

        <div className="container mx-auto px-4 relative z-10">
          <div className="max-w-6xl mx-auto text-center">
            <div className="mb-12">
              <div className="relative w-40 h-40 mx-auto mb-8" style={{ position: 'relative' }}>
                <div className="absolute inset-0 bg-primary/30 rounded-full blur-2xl pulse-glow"></div>
                <Image
                  src="/images/上海生物信息中心表情包v2.png"
                  alt="Join SBC"
                  fill
                  className="object-contain relative z-10"
                />
              </div>
            </div>

            <div className="space-y-8">
              <h2 className="text-5xl lg:text-6xl font-bold">
                <span className="text-foreground/90">JOIN OUR</span>
                <span className="text-primary neon-glow block">TECH COMMUNITY</span>
              </h2>

              <p className="text-xl text-foreground/70 max-w-3xl mx-auto leading-relaxed">
                如果你对生物信息学充满热情，希望与志同道合的伙伴一起交流和学习，
                <span className="text-primary font-semibold">
                  欢迎加入我们的前沿科技学术联盟
                </span>
              </p>

              <div className="flex flex-col sm:flex-row gap-6 justify-center items-center">
                <Button asChild size="lg" className="gap-3 px-8 py-4 text-lg
                         bg-primary hover:bg-primary/90 transition-all duration-300
                         hover:shadow-2xl hover:shadow-primary/30 hover:scale-105 transform">
                  <a href="mailto:contact@sbc.org.cn">
                    <Mail className="h-5 w-5" />
                    联系我们
                  </a>
                </Button>
                <Button asChild variant="outline" size="lg" className="gap-3 px-8 py-4 text-lg
                         border-accent/20 hover:border-accent/40 hover:bg-accent/10
                         hover:shadow-2xl hover:shadow-accent/20 hover:scale-105 transform">
                  <a href="/meetings">
                    <Calendar className="h-5 w-5" />
                    参加活动
                  </a>
                </Button>
              </div>

              <div className="flex flex-wrap gap-4 justify-center items-center">
                {[
                  { icon: "📧", text: "contact@sbc.org.cn" },
                  { icon: "📍", text: "上海生物信息学中心" },
                  { icon: "🎓", text: "欢迎所有年级同学" },
                  { icon: "🚀", text: "一起探索科技未来" }
                ].map((item, index) => (
                  <span
                    key={index}
                    className="px-4 py-2 rounded-full bg-background/80 backdrop-blur-sm
                           border border-primary/20 text-foreground/80 font-mono text-sm
                           hover:bg-primary/10 hover:border-primary/40 transition-all duration-300
                           flex items-center gap-2"
                  >
                    <span>{item.icon}</span>
                    {item.text}
                  </span>
                ))}
              </div>
            </div>

            {/* Tech Decorative Line */}
            <div className="mt-16 flex justify-center">
              <div className="w-32 h-1 bg-gradient-to-r from-primary via-accent to-primary rounded-full"></div>
            </div>
          </div>
        </div>
      </section>

      <Footer />
    </div>
  )
}
